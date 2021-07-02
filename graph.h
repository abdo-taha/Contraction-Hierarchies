#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <limits>
#include <queue>
#include <iostream>
#include <memory>
#include <cassert>
#include <map>
#include <set>
#include <ctime>
#include <stack>

struct path
{
    int cost;
    std::vector<int> nodes;
};



class Graph
{
    typedef int Distance;
    typedef int Vertex;
    typedef std::priority_queue<std::pair<int, int>, std::vector<std::pair<int,int>>, std::greater<std::pair<int, int>>> pr_q;

    // Number of nodes
    int N;
    // Source and target
    int s, t;
    // Estimate of the distance from s to t
    int estimate = INFINITY;
    int mid;
    // Lists of edges outgoing from each node
    std::vector<std::vector<std::pair<int, int>>> outgoing_edges;
    // Lists of edges incoming to each node
    std::vector<std::vector<std::pair<int, int>>> incoming_edges;

    std::vector<std::vector<int>>is_shortcut_out,is_shortcut_in;

    int INFINITY = std::numeric_limits<int>::max() / 4;
    // Levels of nodes for node ordering
    // Ranks of nodes - positions in the node ordering
    std::vector<int> rank,neighbours,level;

    // Distance to node v, bidistance[0][v] - from source in the forward search, bidistance[1][v] - from target
    // in the backward search.
    std::vector<std::vector<Distance>> bidistance;
    std::vector<std::vector<int>> last;
    int m1 = 1,m2 = 2,m3 = 3,m4 = 4;
    
public:
    Graph(int n) {
        // read_stdin();
        set_n(n);
        bidistance.resize(2, std::vector<int>(N, INFINITY));
        last.resize(2,std::vector<int>(n,-1));
    }

    int get_n() { return N;}

    std::vector<std::pair<int, int>>& get_adjacent(int v, bool forward = true) {
        if (forward) {
            return outgoing_edges[v];
        } else {
            return incoming_edges[v];
        }
    }

    void preprocess() {
        pr_q queue;
        for(int i = 0; i < N; ++i){
            if(outgoing_edges[i].size() || incoming_edges[i].size()){
                queue.push(std::make_pair(get_importance(i),i));
            }
        }

        while(queue.size()){
            std::pair<int,int> to_contract = next_v(queue); /// ( importance , v  )
            update_rank(to_contract.second);
            std::vector<Shortcut> SC = shortcuts(to_contract.second);
            do_shortcut(to_contract.second,SC);
            update_neighbour(to_contract.second);
            update_level(to_contract.second);
        }
        dist.clear();
        done.clear();

    }


    path query(int u, int w) {
        estimate = INFINITY;
        s = u;
        t = w;
        pr_q pq[2];
        pq[1].push({0,u});
        pq[0].push({0,w});

        bidistance[1][u] = 0;
        bidistance[0][w] = 0;

        query_id[0][w] = q_id;
        query_id[1][u] = q_id;

        while(pq[0].empty()==false || pq[1].empty()==false){
            if(!pq[1].empty()){
                std::pair<int,int> v = pq[1].top();
                pq[1].pop();
                if(v.first <= estimate){
                    for(std::pair<int,int> &p : outgoing_edges[v.second]) if(rank[v.second] < rank[p.first]){
                        if( query_id[1][p.first] !=  q_id || bidistance[1][p.first] > v.first + p.second ){
                            bidistance[1][p.first] = v.first + p.second;
                            pq[1].push({v.first + p.second,p.first});
                            last[1][p.first] = v.second;
                            query_id[1][p.first] = q_id;
                        }
                    }
                }
                if(query_id[0][v.second] == q_id && bidistance[0][v.second]+bidistance[1][v.second] <estimate ){
                    estimate = bidistance[0][v.second]+bidistance[1][v.second];
                    mid = v.second;
                }
            }
            if(!pq[0].empty()){
                std::pair<int,int> v = pq[0].top();
                pq[0].pop();
                if(v.first <= estimate){
                    for(std::pair<int,int> &p : incoming_edges[v.second])if(rank[v.second] < rank[p.first]){
                        if( query_id[0][p.first] != q_id   || bidistance[0][p.first] > v.first + p.second ){
                            bidistance[0][p.first] = v.first + p.second;
                            pq[0].push({v.first + p.second,p.first});
                            last[0][p.first] = v.second;
                            query_id[0][p.first] = q_id;
                        }
                    }
                }
                if(query_id[1][v.second] == q_id && bidistance[0][v.second]+bidistance[1][v.second] < estimate) {
                    estimate  = bidistance[0][v.second]+bidistance[1][v.second];
                    mid = v.second;
                }
            }
        }

        ++q_id;
        path ans;
        ans.cost = -1;
        if(estimate < INFINITY) {
            ans.cost = estimate;
            ans.nodes = get_path(mid);
        }
        return ans;
    }

private:
    /// -------------------------- preprocessing -------------------------------///
    struct Shortcut {
        int from;
        int to;
        int cost;
        Shortcut(int from, int to, int cost):from(from),to(to),cost(cost){}
    };

    int get_importance (int v){
        return m1 * ( outgoing_edges[v].size() * incoming_edges[v].size() - outgoing_edges[v].size() - incoming_edges[v].size())
        +   m2 *(outgoing_edges[v].size() + incoming_edges[v].size())
        + m3 * neighbours[v]
        + m4* level[v];
    }

    std::pair<int,int> next_v(pr_q &queue){
        std::pair<int,int> to_contract = queue.top();
        while(true){
            queue.pop();
            to_contract.first = get_importance(to_contract.second);

            if(queue.empty() ||to_contract.first <= queue.top().first){
                break;
            }
            queue.push(to_contract);
            to_contract = queue.top();
        }
        return to_contract;
    }

    int rnk = 1;
    void update_rank(int v){
        rank[v] = rnk++;
    }

    void update_level(int v){
        for(std::pair<int,int> &before : incoming_edges[v]) level[before.first] = std::max(level[before.first],level[v]+1);
        for(std::pair<int,int> &nxt : outgoing_edges[v]) level[nxt.first] = std::max(level[nxt.first],level[v]+1);
    }

    void update_neighbour(int v){
        for(std::pair<int,int> &before : incoming_edges[v]) ++neighbours[before.first];
        for(std::pair<int,int> &nxt : outgoing_edges[v]) ++neighbours[nxt.first];
    }

    std::vector<Shortcut> shortcuts(int v){
        std::vector<Shortcut> shorts;
        int mxb = 0 , mxn = 0;
        for(std::pair<int,int> &before : incoming_edges[v]) if(rank[before.first] == INFINITY) mxb = std::max(mxb,before.second);
        for(std::pair<int,int> &nxt : outgoing_edges[v]) if(rank[nxt.first] == INFINITY) mxn = std::max(mxn,nxt.second);
        int mx = mxb + mxn;

        for(std::pair<int,int> &before : incoming_edges[v]){
            if(rank[before.first]  != INFINITY) continue;
            dijkstra(before.first,v,mx,5000);
            for(std::pair<int,int> &nxt : outgoing_edges[v])if(before.first!=nxt.first){
                if(rank[nxt.first] != INFINITY) continue;

                if(dist[nxt.first] > before.second+nxt.second){
                    shorts.push_back( Shortcut(before.first,nxt.first,before.second+nxt.second));
                }
            }
            clear_done();
        }
        return shorts;
    }

    std::vector<int> dist; // for preprocessing (check witness path)
    std::set<int> done; // used in last witness

    // witness search
    void dijkstra(int from , int v, int mx, int hops){ 
        using T =  std::pair<std::pair<int,int>,int>; /// ( (dist ,hops) , v)
        std::priority_queue<T,std::vector<T>,std::greater<T>> pq;
        dist[from] = 0;
        done.insert(from);
        pq.push(std::make_pair(std::make_pair(0,0),from));
        while(pq.size()){
            int cur = pq.top().second, cost = pq.top().first.first,hps =pq.top().first.second ;
            pq.pop();
            if(cost != dist[cur]) continue;
            if(cost >= mx ) return ; /// ||hps> hops
            for(std::pair<int,int> &nxt : outgoing_edges[cur]){
                if( rank[nxt.first] == INFINITY && nxt.first != v &&  dist[nxt.first] > cost + nxt.second ){
                    dist[nxt.first] = cost+nxt.second;
                    pq.push( std::make_pair(std::make_pair(cost+nxt.second,hps+1),nxt.first) );
                    done.insert(nxt.first);
                }
            }
        }
    }

    void clear_done(){
        for(int x : done) dist[x] = INFINITY;
        done.clear();
    }

    void do_shortcut(int v, std::vector<Shortcut>& shortcuts) {
        for(Shortcut &sc : shortcuts){
            add_directed_edge(sc.from,sc.to,sc.cost,v);
        }
    }

    ///  ---------------------------------- queries ------------------------------------------------------///


    std::vector<int> query_id[2];
    int q_id = 1;
    std::vector<int> get_path(int v){
        std::cout << s << " " << v << " " << t << "\n";
        std::stack<int> st;
        int back = v;
        while (back != s)
        {
            int before = last[1][back];
            std::vector<int> tmp = unshort(before,back,outgoing_edges,is_shortcut_out);
            for (int i = tmp.size()-1; i > 0; i--)
            {
                st.push(tmp[i]);
            }
            back = before;
        }
        std::vector<int> ans;
        ans.push_back(s);
        while (st.size())
        {
            ans.push_back(st.top());
            st.pop();
        }
        

        while (v != t)
        {
            int before = last[0][v];
            std::vector<int> tmp = unshort(before,v,incoming_edges,is_shortcut_in);
            for (int i = tmp.size()-2; i >= 0; i--)
            {
                ans.push_back(tmp[i]);
            }
            v = before;
        }
        return ans;
        
    }

    std::vector<int> unshort(int from ,int to,std::vector<std::vector<std::pair<int, int>>> &edges,std::vector<std::vector<int>> &shortcut){
        std::vector<int> a,b;
        int md = get_shortcut(edges[from],shortcut[from],to);
        if( md == -1){
            a.push_back(from);
            a.push_back(to);
            return a;
        }else{
            a = unshort(from,md,edges,shortcut);
            b = unshort(md,to,edges,shortcut);
            for (int i = 1; i < b.size(); i++)
            {
                a.push_back(b[i]);
            }
            return a;
             
        }

    }
    
    /// -------------------------------------  input --------------------------------------------///
    void set_n(int n) {
        N = n;
        outgoing_edges.resize(n);
        incoming_edges.resize(n);
        is_shortcut_in.resize(n);
        is_shortcut_out.resize(n);
        rank.assign(n,INFINITY);
        level.assign(n,0);
        neighbours.assign(n,0);
        dist.assign(n,INFINITY);
        query_id[0].assign(n,0);
        query_id[1].assign(n,0);
    }

    void add_edge_to_list(std::vector<std::pair<int,int>>& list,std::vector<int>&is_shortcut, int w, int c,int shortcut){
        for (int i = 0; i < list.size(); ++i) {
            std::pair<int, int>& p = list[i];
            if (p.first == w) {
                if (p.second > c) {
                    p.second = c;
                    is_shortcut[i] = shortcut;
                }
                return;
            }
        }
        list.push_back(std::make_pair(w, c));
        is_shortcut.push_back(shortcut);
    }
    
    int get_shortcut(std::vector<std::pair<int,int>>& list,std::vector<int>&is_shortcut,int to){
        for (int i = 0; i < list.size(); ++i) {
            std::pair<int, int>& p = list[i];
            if (p.first == to) {
                    return is_shortcut[i];
            }
        }
        return -2;/// .....
    }

    void add_directed_edge(int u, int v, int c,int shortcut) {
        add_edge_to_list(outgoing_edges[u],is_shortcut_out[u], v, c,shortcut);
        add_edge_to_list(incoming_edges[v],is_shortcut_in[v],u, c,shortcut);
    }

public :
    void add_edge(int from, int to, int cost) {
        add_directed_edge(from, to, cost,-1);
    }
    void change_importance(int edge_diff,int short_cover,int contracted_neigbor,int nlevel){
        m1 = edge_diff;
        m2 = short_cover;
        m3 = contracted_neigbor;
        m4 = nlevel;
    }
};
