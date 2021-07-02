# Contraction Hierarchies
advanced algorithm for comupting shortest path

The code was tested on a [road network graph of New York](http://www.dis.uniroma1.it/challenge9/download.shtml) (has about 250K nodes, 700K edges) and takes on average 31s for preprocessing and average query time of 2.5ms on my machine.

[wikipidia](https://en.wikipedia.org/wiki/Contraction_hierarchies)
[source](https://d3c33hcgiwev3.cloudfront.net/_4b3e617775b52e3c72a457a310a2be43_19_advanced_shortest_paths_3_contraction_hierarchies.pdf?Expires=1625270400&Signature=RMaGLTZ950PCY8x4joLlIMetizXeQVpgA8Du4pPLcFmnrXoDuH2RmA7MPhFNmResVOBqy0TYWsY~Tn-W-Bcv-ufN68IZVxQ92~SVyGFYSFQZPWo06cogMRqbEZYvxLrZbl1lvhdhrVFxrtcxwmCggSTWOUs2~ha~nMejUNYzYqg_&Key-Pair-Id=APKAJLTNE6QMUY6HBC5A)

## usage 
-download header file.
-include it.

-use class graph(number of nodes).

-use add_edge(from,to,cost) nodes should be zero based.

-preprocess().

-query(from,to).

-query return struct path which have cost(int) nodes(vector<int>).

-- change_importance can be used before preprocessing to change importance constants.


```cpp
#include "graph.h"
int main() {
    Graph g(100);
    for (int i = 0; i < 100; i++)
    {
        g.add_edge(from-1,to-1,cost);
    }
    g.preprocess();
    path tmp = g.query(from-1,to-1);
    printf("%d\n", tmp.cost);
    for(int x : tmp.nodes) printf("%d ",x+1); 
    return 0;
}
```
# important

use -O flag for compiler to optimize code