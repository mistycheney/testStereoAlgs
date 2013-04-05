#ifndef __MAXFLOW_H__
#define __MAXFLOW_H__


typedef long flowtype;

struct Edge
{
	long		from;
	long		to;
	long		cap;
};

#define BLOCK_SIZE 512

struct EdgeList
{
	struct Edge edge[BLOCK_SIZE];
	struct EdgeList *next;
	int    num;
};

struct Graph
{
	long		node_min, node_max;
	long		m;
	long		source;
	long		sink;
	struct EdgeList	*first;
	struct EdgeList *current;
};

struct Graph * init_graph(long source, long sink);
void add_edge(struct Graph *graph, long from, long to, long cap);
flowtype maxflow(struct Graph *graph, int *cut);

#endif
