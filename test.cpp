#include "graph.h"

int main() {

    std::ios::sync_with_stdio(0);
    auto c1 = clock();

    freopen("newyork.txt","r",stdin);
    freopen("o.txt","w",stdout);
    int size,edges;
    scanf("%d %d" , &size, &edges);
    Graph g(size);
    for (int i = 0; i < edges; i++)
    {
        int from,to,cost;
        scanf("%d %d %d" ,&from,&to,&cost);
        g.add_edge(from-1,to-1,cost);
    }
    
    g.preprocess();
    auto c2 = clock();
    int t;
    double t_q = 0;
    assert(scanf("%d", &t) == 1);
    for (int i = 0; i < t; ++i) {
        int u, v;
        assert(scanf("%d %d", &u, &v) == 2);
        auto c3 = clock();
        path tmp = g.query(u-1,v-1);
        auto c4 = clock();
        printf("%d\n", tmp.cost);
        for(int x : tmp.nodes) printf("%d ",x+1);
        printf("\n");
        t_q += 1.0 * (c4-c3);
    }
    auto c3 = clock();
    printf("pre : %f    q: %f\n",(double) (1.0*(c2-c1) / CLOCKS_PER_SEC) , (double)(t_q/ CLOCKS_PER_SEC / t));
}
