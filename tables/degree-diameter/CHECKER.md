# Checker: graphs and digraphs of given degree and diameter

Open for the vertex-symmetric digraphs table: Cayley or coset digraphs of PSL(2, q) or PGL(2, q), times a cyclic group Z_m if wanted, given by a JSON spec (the group, the generators as 2 x 2 matrices over GF(q) with an element of Z_m, the claimed out-degree d, diameter k and order n). A Cayley or coset digraph is vertex-transitive by construction, so its diameter is the eccentricity of one vertex. The checkers are the two programs of github.com/arhancanli/quantavane-vsd, folder `checkers/`, written separately:

- `vsd_check.py`: rebuilds the group as permutations of the projective line times Z_m by its own search, builds the digraph from the definition, checks degrees, order and the base vertex's eccentricity, and for n <= 20,000 also the eccentricity of every vertex by reachability bitmasks;
- `vsd_check.c`: reads only the explicit arc and label files and checks them directly (`vsd_check_big.py` for the two largest orders).

```sh
python3 -I checkers/vsd_check.py SPEC.json [ARCS LABELS]
```

The other three tables here (Cayley graphs, bipartite graphs, arc-transitive graphs) are listed but not open for submissions yet: a checker for them must also verify the symmetry, which these programs do not do for general graphs.

## Runs on our own certificates (19 November 2026)

- vsd_check.py on the public specs of the (2,10) digraph on 660 vertices and the (3,10) digraph on 12,180 vertices: OK (out-degree, order and diameter as claimed).
- Each packet in quantavane-vsd keeps the output of both programs and of red controls (checks.txt).
