#include "maxflow.h"

/* Comment out the next line if you are getting compile problems (see README-maxflow.txt): */
// #define PRF_DOWNLOADED

#ifndef PRF_DOWNLOADED

/*  These stubs are defined if the push-relabel flow code has not been downloaded
    (see README.txt) */

struct Graph * init_graph(long source, long sink)
{
    return (struct Graph *) 0;
}

void add_edge(struct Graph *graph, long from, long to, long cap) {}

flowtype maxflow (struct Graph *graph, int *cut)
{
    return 0;
}

#else

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXLONG ((1 << 30) - 1)

#if (!defined(H_PRF) || !defined(M_PRF) || !defined(Q_PRF) || !defined(F_PRF) || !defined(DF))
#define H_PRF
#endif
/*
    Define one of the values here:
    H_PRF	push-relabel algorithm (highest level). 
    M_PRF	push-relabel algorithm (highest level), no gaps. 
    Q_PRF	push-relabel algorithm (queue). 
    F_PRF	push-relabel algorithm (queue), no gaps. 
    DF		Dinitz algorithm. 
*/

/* statistic variables */
long n_push  = 0;         /* number of pushes */
long n_rel   = 0;         /* number of relabels */
long n_up    = 0;         /* number of updates */
long n_gap   = 0;         /* number of gaps */
long n_gnode = 0;         /* number of nodes after gap */
float t, t2;

float timer () { return 0.0; }

#ifdef H_PRF
/* Comment out #define PRF_DOWNLOADED line at the top if you are getting compile problems: */
#include "types_pr.h"
#include "parser_flg.c"
#include "h_prf.c"
#endif

#ifdef M_PRF
#include "types_pr.h"
#include "parser_flg.c"
#include "m_prf.c"
#endif

#ifdef Q_PRF
#include "types_qpr.h"
#include "parser_flg.c"
#include "q_prf.c"
#endif

#ifdef F_PRF
#include "types_fpr.h"
#include "parser_flg.c"
#include "f_prf.c"
#endif

#ifdef DF
#include "types_df.h"
#include "parser_flg.c"
#include "df.c"
#endif

void mem_error()
{
	fprintf ( stderr, "Allocating error\n");
	exit ( 1 );
}

flowtype maxflow (struct Graph *graph, int *cut)
{
	arc *arp;
	long *cap;
	node *ndp, *source, *sink, *i;
	long n, m, nmin;
	long ni;
	int  cc;
	flowtype flow = 0;

	cc = parse( graph, &n, &m, &ndp, &arp, &cap, &source, &sink, &nmin );
	if ( cc ) mem_error();
#ifdef DF
	cc = dflow ( n, ndp, arp, source, sink, &flow );
#else
	cc = prflow ( n, ndp, arp, cap, source, sink, &flow );
#endif

    /* The following 3 (ifdefed) lines copied from
       Vladimir Kolmogorov's version of h_prf.c: */
#ifndef DF
    def_ranks();
#endif
    free ( queue );
#if defined(H_PRF) || defined (M_PRF) || defined(Q_PRF)
    free ( layers );
#endif

	if ( cc ) mem_error();

	i = ndp;
	for ( ni = graph->node_min; ni < graph->node_max; ni++ )
	{
		cut[ni] = ((i++)->rank == source->rank) ? 0 : 1;
	}

	free(graph);
	free(ndp - nmin);
	free(arp);
	free(cap);

	return flow;
}

struct Graph * init_graph(long source, long sink)
{
	struct Graph * graph = (struct Graph *) malloc(sizeof(struct Graph));
	if(graph == NULL) mem_error();

	graph->m = 0;
	graph->node_min = (source < sink) ? source : sink;
	graph->node_max = (source > sink) ? source+1 : sink+1;
	graph->source = source; graph->sink = sink;
	graph->first = (struct EdgeList *) malloc(sizeof(struct EdgeList));
	if(graph->first == NULL) mem_error();
	graph->first->next = NULL;
	graph->first->num = 0;
	graph->current = graph->first;

	return graph;
}

void add_edge(struct Graph *graph, long from, long to, long cap)
{
	if(graph->current->num == BLOCK_SIZE)
	{
		graph->current->next = (struct EdgeList *) malloc(sizeof(struct EdgeList));
		if(graph->current->next == NULL) mem_error();
		graph->current = graph->current->next;
		graph->current->next = NULL;
		graph->current->num = 0;
	}
	graph->m++;
	if(graph->node_min > from) graph->node_min = from;
	if(graph->node_max < from+1) graph->node_max = from+1;
	if(graph->node_min > to) graph->node_min = to;
	if(graph->node_max < to+1) graph->node_max = to+1;
	graph->current->edge[graph->current->num].from = from;
	graph->current->edge[graph->current->num].to = to;
	graph->current->edge[graph->current->num].cap = cap;
	graph->current->num++;
}

#endif /* PRF_DOWNLOADED */
