//! Minimal GFA 1.0 parser sufficient for impg-produced subgraphs.
//!
//! Only handles S (segment), L (link), and P (path) lines. Strands on links
//! are assumed all-forward (impg-emitted GFAs from `impg query` produce that).

use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path as FsPath;

pub type NodeId = u64;

/// A path through the graph (one haplotype's traversal).
#[derive(Debug, Clone)]
pub struct Path {
    pub name: String,
    pub nodes: Vec<NodeId>,
}

/// A pangenome subgraph parsed from a GFA file.
#[derive(Debug)]
pub struct Graph {
    /// Sequence per node (one or more bases; impg emits one base per node).
    pub seq: HashMap<NodeId, Vec<u8>>,
    /// Forward adjacency: node → list of immediate successors.
    pub forward: HashMap<NodeId, Vec<NodeId>>,
    /// Backward adjacency: node → list of immediate predecessors.
    pub backward: HashMap<NodeId, Vec<NodeId>>,
    /// Paths in the graph.
    pub paths: Vec<Path>,
}

impl Graph {
    pub fn parse(path: &FsPath) -> Result<Self> {
        let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let r = BufReader::new(f);

        let mut seq: HashMap<NodeId, Vec<u8>> = HashMap::new();
        let mut forward: HashMap<NodeId, Vec<NodeId>> = HashMap::new();
        let mut backward: HashMap<NodeId, Vec<NodeId>> = HashMap::new();
        let mut paths: Vec<Path> = Vec::new();

        for (i, line) in r.lines().enumerate() {
            let line = line.with_context(|| format!("reading line {}", i + 1))?;
            if line.is_empty() {
                continue;
            }
            let mut fields = line.split('\t');
            let tag = fields.next().unwrap_or("");
            match tag {
                "H" => continue,
                "S" => {
                    let id: NodeId = fields
                        .next()
                        .ok_or_else(|| anyhow::anyhow!("S line missing id (line {})", i + 1))?
                        .parse()
                        .with_context(|| format!("parsing S id on line {}", i + 1))?;
                    let s = fields
                        .next()
                        .ok_or_else(|| anyhow::anyhow!("S line missing seq (line {})", i + 1))?;
                    seq.insert(id, s.as_bytes().to_vec());
                }
                "L" => {
                    let from: NodeId = fields
                        .next()
                        .and_then(|f| f.parse().ok())
                        .ok_or_else(|| anyhow::anyhow!("L line bad from on line {}", i + 1))?;
                    let _from_strand = fields.next();
                    let to: NodeId = fields
                        .next()
                        .and_then(|f| f.parse().ok())
                        .ok_or_else(|| anyhow::anyhow!("L line bad to on line {}", i + 1))?;
                    let _to_strand = fields.next();
                    forward.entry(from).or_default().push(to);
                    backward.entry(to).or_default().push(from);
                }
                "P" => {
                    let name = fields
                        .next()
                        .ok_or_else(|| anyhow::anyhow!("P line missing name on line {}", i + 1))?
                        .to_string();
                    let path_str = fields
                        .next()
                        .ok_or_else(|| anyhow::anyhow!("P line missing path on line {}", i + 1))?;
                    let mut nodes = Vec::new();
                    for token in path_str.split(',') {
                        // Strip optional strand suffix (+ or -).
                        let token = token.trim_end_matches(['+', '-']);
                        let id: NodeId = token
                            .parse()
                            .with_context(|| format!("parsing path node on line {}", i + 1))?;
                        nodes.push(id);
                    }
                    paths.push(Path { name, nodes });
                }
                _ => continue, // ignore other tag types (E, W, etc.)
            }
        }

        if paths.is_empty() {
            bail!("no paths found in {}", path.display());
        }
        Ok(Graph { seq, forward, backward, paths })
    }

    /// Successors of a node, empty slice if none.
    pub fn successors(&self, n: NodeId) -> &[NodeId] {
        self.forward.get(&n).map_or(&[], |v| v.as_slice())
    }

    /// Predecessors of a node, empty slice if none.
    pub fn predecessors(&self, n: NodeId) -> &[NodeId] {
        self.backward.get(&n).map_or(&[], |v| v.as_slice())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn parse_gfa(text: &str) -> Result<Graph> {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(text.as_bytes()).unwrap();
        f.flush().unwrap();
        Graph::parse(f.path())
    }

    /// Diamond: 1 → {2, 3} → 4, one path 1,2,4 and another 1,3,4.
    fn diamond_gfa() -> &'static str {
        "H\tVN:Z:1.0\n\
         S\t1\tA\n\
         S\t2\tC\n\
         S\t3\tG\n\
         S\t4\tT\n\
         L\t1\t+\t2\t+\t0M\n\
         L\t1\t+\t3\t+\t0M\n\
         L\t2\t+\t4\t+\t0M\n\
         L\t3\t+\t4\t+\t0M\n\
         P\thapA\t1+,2+,4+\t*\n\
         P\thapB\t1+,3+,4+\t*\n"
    }

    #[test]
    fn parse_diamond_builds_segments_links_paths() {
        let g = parse_gfa(diamond_gfa()).unwrap();
        assert_eq!(g.seq.len(), 4);
        assert_eq!(g.seq[&1], b"A");
        assert_eq!(g.seq[&4], b"T");
        assert_eq!(g.paths.len(), 2);
        assert_eq!(g.paths[0].name, "hapA");
        assert_eq!(g.paths[0].nodes, vec![1, 2, 4]);
        assert_eq!(g.paths[1].nodes, vec![1, 3, 4]);
    }

    #[test]
    fn successors_predecessors_match_links() {
        let g = parse_gfa(diamond_gfa()).unwrap();
        let mut succ_of_1: Vec<_> = g.successors(1).to_vec();
        succ_of_1.sort_unstable();
        assert_eq!(succ_of_1, vec![2, 3]);
        let mut pred_of_4: Vec<_> = g.predecessors(4).to_vec();
        pred_of_4.sort_unstable();
        assert_eq!(pred_of_4, vec![2, 3]);
        // Terminal nodes have no successors / no predecessors.
        assert!(g.successors(4).is_empty());
        assert!(g.predecessors(1).is_empty());
    }

    #[test]
    fn missing_paths_errors() {
        let text = "S\t1\tA\nS\t2\tC\nL\t1\t+\t2\t+\t0M\n";
        let err = parse_gfa(text).unwrap_err();
        assert!(format!("{err}").contains("no paths found"));
    }

    #[test]
    fn empty_lines_and_unknown_tags_are_ignored() {
        let text = "H\tVN:Z:1.0\n\
                    \n\
                    S\t1\tA\n\
                    S\t2\tT\n\
                    E\t1\t+\t2\t+\n\
                    W\thapA\t1\thapA-coords\n\
                    L\t1\t+\t2\t+\t0M\n\
                    P\thapA\t1+,2+\t*\n";
        let g = parse_gfa(text).unwrap();
        assert_eq!(g.seq.len(), 2);
        assert_eq!(g.paths.len(), 1);
    }

    #[test]
    fn path_strand_suffix_is_stripped() {
        // Path tokens with mixed strand markers must parse as plain ids.
        let text = "S\t1\tA\nS\t2\tT\nL\t1\t+\t2\t+\t0M\nP\thapA\t1-,2+\t*\n";
        let g = parse_gfa(text).unwrap();
        assert_eq!(g.paths[0].nodes, vec![1, 2]);
    }

    #[test]
    fn successors_predecessors_missing_node_returns_empty() {
        let g = parse_gfa(diamond_gfa()).unwrap();
        // 999 doesn't exist in the graph; lookup must not panic.
        assert!(g.successors(999).is_empty());
        assert!(g.predecessors(999).is_empty());
    }

    #[test]
    fn bad_segment_id_errors() {
        let text = "S\tnotanid\tA\nP\tp\t1+\t*\n";
        assert!(parse_gfa(text).is_err());
    }
}
