//! Bubble enumeration on a pangenome graph.
//!
//! A bubble here is a pair (source, sink) such that every successor of `source`
//! reaches `sink` (forward) before reconverging elsewhere. The "branches"
//! are the internal node sequences between source and sink, exclusive of both.
//!
//! v0.1: top-level bubbles only — bubbles nested inside a branch are not
//! enumerated separately. Sufficient for impg-emitted GFAs of short regions.

use std::collections::hash_map::Entry;
use std::collections::HashMap;

use crate::gfa::{Graph, NodeId};

/// A bubble found in the graph.
#[derive(Debug, Clone)]
pub struct Bubble {
    pub source: NodeId,
    pub sink: NodeId,
    /// One internal-node list per branch (in order of `Graph::successors(source)`).
    /// An empty branch means the branch goes directly from source to sink.
    pub branches: Vec<Vec<NodeId>>,
}

impl Bubble {
    pub fn n_branches(&self) -> usize {
        self.branches.len()
    }
}

/// Find the bubble rooted at `source` if one closes within `max_depth` steps.
///
/// Uses parallel BFS from each immediate successor, tracking reachable sets
/// per branch. The sink is the node reachable from every branch with the
/// smallest maximum depth across branches.
pub fn find_bubble(graph: &Graph, source: NodeId, max_depth: usize) -> Option<Bubble> {
    let succs = graph.successors(source);
    if succs.len() < 2 {
        return None;
    }
    let n = succs.len();

    let mut depth: Vec<HashMap<NodeId, usize>> = vec![HashMap::new(); n];
    let mut parents: Vec<HashMap<NodeId, NodeId>> = vec![HashMap::new(); n];
    let mut frontier: Vec<Vec<NodeId>> = vec![Vec::new(); n];

    for (i, &s) in succs.iter().enumerate() {
        depth[i].insert(s, 0);
        frontier[i].push(s);
    }

    for step in 0..max_depth {
        // Candidates: nodes present in every branch's reachable set.
        let candidates: Vec<NodeId> = depth[0]
            .keys()
            .copied()
            .filter(|n_| depth[1..].iter().all(|m| m.contains_key(n_)))
            .collect();
        if !candidates.is_empty() {
            // Pick the sink as the node minimizing max-depth across branches
            // (= earliest reachable by the slowest branch).
            let sink = *candidates
                .iter()
                .min_by_key(|n_| depth.iter().map(|d| d[n_]).max().unwrap_or(usize::MAX))
                .unwrap();
            let branches = reconstruct_branches(&parents, succs, sink);
            return Some(Bubble { source, sink, branches });
        }

        // Expand all frontiers by one step.
        let cur_depth = step + 1;
        let mut next_frontier: Vec<Vec<NodeId>> = vec![Vec::new(); n];
        let mut any_progress = false;
        for i in 0..n {
            for &node in &frontier[i] {
                for &m in graph.successors(node) {
                    if let Entry::Vacant(e) = depth[i].entry(m) {
                        e.insert(cur_depth);
                        parents[i].insert(m, node);
                        next_frontier[i].push(m);
                        any_progress = true;
                    }
                }
            }
        }
        frontier = next_frontier;
        if !any_progress {
            return None;
        }
    }
    None
}

fn reconstruct_branches(
    parents: &[HashMap<NodeId, NodeId>],
    succs: &[NodeId],
    sink: NodeId,
) -> Vec<Vec<NodeId>> {
    let n = succs.len();
    let mut branches = Vec::with_capacity(n);
    for i in 0..n {
        // If the immediate successor IS the sink, branch is empty.
        if succs[i] == sink {
            branches.push(Vec::new());
            continue;
        }
        // Walk back from sink through parents[i] until we hit succs[i].
        let mut path = vec![sink];
        let mut cur = sink;
        while let Some(&p) = parents[i].get(&cur) {
            cur = p;
            path.push(cur);
            if cur == succs[i] {
                break;
            }
        }
        path.reverse();
        // path now starts at succs[i] and ends at sink. Drop the sink so the
        // branch is just the internal nodes.
        if path.last() == Some(&sink) {
            path.pop();
        }
        // Sanity: path[0] should be succs[i]. If not, the BFS state was
        // inconsistent — fall back to an empty branch rather than panic.
        if path.first() != Some(&succs[i]) {
            branches.push(Vec::new());
            continue;
        }
        branches.push(path);
    }
    branches
}

/// All nodes with ≥2 outgoing edges, sorted by node id for determinism.
/// These are the candidate sources for classification — every such node
/// represents a bubble in the broad sense (a branching point in the graph),
/// whether or not its BFS sink converges.
pub fn enumerate_sources(graph: &Graph) -> Vec<NodeId> {
    let mut sources: Vec<NodeId> = graph
        .forward
        .iter()
        .filter(|(_, v)| v.len() >= 2)
        .map(|(k, _)| *k)
        .collect();
    sources.sort_unstable();
    sources
}

/// Enumerate sources, returning a Bubble per source. When BFS converges,
/// the bubble has a real sink and reconstructed branches. When it doesn't
/// (e.g. interior structure that doesn't close within `max_depth`), we still
/// emit a bubble with `sink = source` and empty `branches` so the classifier
/// can still walk from the source and decide a type (typically Complex).
pub fn enumerate_bubbles(graph: &Graph, max_depth: usize) -> Vec<Bubble> {
    enumerate_sources(graph)
        .into_iter()
        .map(|s| {
            find_bubble(graph, s, max_depth).unwrap_or(Bubble {
                source: s,
                sink: s,
                branches: Vec::new(),
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gfa::Path;
    use std::collections::HashMap;

    /// Build a Graph directly from forward edges + path names.
    /// Each segment has a single byte sequence equal to its id mod 26 + 'A'.
    fn graph(forward_edges: &[(NodeId, NodeId)], paths: &[(&str, &[NodeId])]) -> Graph {
        let mut seq: HashMap<NodeId, Vec<u8>> = HashMap::new();
        let mut forward: HashMap<NodeId, Vec<NodeId>> = HashMap::new();
        let mut backward: HashMap<NodeId, Vec<NodeId>> = HashMap::new();
        for &(a, b) in forward_edges {
            forward.entry(a).or_default().push(b);
            backward.entry(b).or_default().push(a);
            seq.entry(a).or_insert_with(|| vec![b'A' + (a % 26) as u8]);
            seq.entry(b).or_insert_with(|| vec![b'A' + (b % 26) as u8]);
        }
        let paths = paths
            .iter()
            .map(|(n, ns)| Path { name: (*n).to_string(), nodes: ns.to_vec() })
            .collect();
        Graph { seq, forward, backward, paths }
    }

    #[test]
    fn enumerate_sources_picks_only_branching_nodes() {
        // 1 → {2, 3}, 2 → 4, 3 → 4 — only node 1 has ≥2 outgoing.
        let g = graph(
            &[(1, 2), (1, 3), (2, 4), (3, 4)],
            &[("a", &[1, 2, 4]), ("b", &[1, 3, 4])],
        );
        assert_eq!(enumerate_sources(&g), vec![1]);
    }

    #[test]
    fn enumerate_sources_returns_sorted_deterministic() {
        // Two branch points: 1 → {2, 3} and 5 → {6, 7}.
        let g = graph(
            &[(1, 2), (1, 3), (5, 6), (5, 7), (2, 5), (3, 5), (6, 8), (7, 8)],
            &[("a", &[1, 2, 5, 6, 8]), ("b", &[1, 3, 5, 7, 8])],
        );
        assert_eq!(enumerate_sources(&g), vec![1, 5]);
    }

    #[test]
    fn find_bubble_diamond() {
        // 1 → {2, 3} → 4.
        let g = graph(
            &[(1, 2), (1, 3), (2, 4), (3, 4)],
            &[("a", &[1, 2, 4]), ("b", &[1, 3, 4])],
        );
        let b = find_bubble(&g, 1, 8).expect("diamond should yield a bubble");
        assert_eq!(b.source, 1);
        assert_eq!(b.sink, 4);
        assert_eq!(b.n_branches(), 2);
        // Each branch is exactly one internal node (2 or 3).
        let mut internals: Vec<NodeId> = b.branches.iter().flatten().copied().collect();
        internals.sort_unstable();
        assert_eq!(internals, vec![2, 3]);
    }

    #[test]
    fn find_bubble_returns_none_for_non_branching_source() {
        let g = graph(&[(1, 2), (2, 3)], &[("a", &[1, 2, 3])]);
        assert!(find_bubble(&g, 1, 8).is_none());
    }

    #[test]
    fn find_bubble_returns_none_when_branches_never_reconverge() {
        // 1 → {2, 3}; 2 ends, 3 ends. No common descendant.
        let g = graph(&[(1, 2), (1, 3)], &[("a", &[1, 2]), ("b", &[1, 3])]);
        assert!(find_bubble(&g, 1, 8).is_none());
    }

    #[test]
    fn find_bubble_respects_max_depth() {
        // Long bubble: 1 → {2, 3}; chains 2 → 5 → 7 and 3 → 6 → 7. Sink at depth 2.
        // max_depth is the iteration cap; the candidate check at iteration N
        // sees frontier depth N, so a depth-2 sink needs max_depth >= 3.
        let g = graph(
            &[(1, 2), (1, 3), (2, 5), (3, 6), (5, 7), (6, 7)],
            &[("a", &[1, 2, 5, 7]), ("b", &[1, 3, 6, 7])],
        );
        // max_depth=2 cannot reach the depth-2 candidate check → no bubble.
        assert!(find_bubble(&g, 1, 2).is_none());
        // max_depth=3 includes the iteration where the sink becomes shared.
        let b = find_bubble(&g, 1, 3).expect("bubble closes by depth 2");
        assert_eq!(b.sink, 7);
    }

    #[test]
    fn find_bubble_handles_empty_branch_to_sink() {
        // 1 → {2, 4}, 2 → 4 — branch via 2 has internal node 2, branch via 4
        // is empty (direct edge to sink).
        let g = graph(&[(1, 2), (1, 4), (2, 4)], &[("a", &[1, 2, 4]), ("b", &[1, 4])]);
        let b = find_bubble(&g, 1, 8).expect("bubble exists");
        assert_eq!(b.sink, 4);
        // Exactly one branch is empty (direct), the other has [2].
        let empty_count = b.branches.iter().filter(|br| br.is_empty()).count();
        let single_node_count = b.branches.iter().filter(|br| br.len() == 1).count();
        assert_eq!(empty_count, 1);
        assert_eq!(single_node_count, 1);
    }

    #[test]
    fn enumerate_bubbles_emits_one_per_source_with_fallback() {
        // Two sources: diamond at 1 (converges), open at 5 (no reconvergence).
        let g = graph(
            &[(1, 2), (1, 3), (2, 4), (3, 4), (5, 6), (5, 7)],
            &[("a", &[1, 2, 4, 5, 6]), ("b", &[1, 3, 4, 5, 7])],
        );
        let bubbles = enumerate_bubbles(&g, 8);
        assert_eq!(bubbles.len(), 2);
        // Diamond at source 1.
        let b1 = bubbles.iter().find(|b| b.source == 1).unwrap();
        assert_eq!(b1.sink, 4);
        assert_eq!(b1.n_branches(), 2);
        // Open bubble at source 5 → fallback (sink == source, no branches).
        let b5 = bubbles.iter().find(|b| b.source == 5).unwrap();
        assert_eq!(b5.sink, 5);
        assert_eq!(b5.n_branches(), 0);
    }

    #[test]
    fn three_branch_bubble() {
        // 1 → {2, 3, 4} → 5.
        let g = graph(
            &[(1, 2), (1, 3), (1, 4), (2, 5), (3, 5), (4, 5)],
            &[("a", &[1, 2, 5]), ("b", &[1, 3, 5]), ("c", &[1, 4, 5])],
        );
        let b = find_bubble(&g, 1, 4).expect("three-way bubble");
        assert_eq!(b.sink, 5);
        assert_eq!(b.n_branches(), 3);
    }
}
