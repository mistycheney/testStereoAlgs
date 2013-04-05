/*
  parse (...) :
       1. Reads maximal flow problem in extended DIMACS format.
       2. Prepares internal data representation.
       
   types: 'arc' and 'node' must be predefined

   type arc  must contain fields 'head', 'sister', 'next', 'r_cap'

   typedef 
     struct arc_st
       {
         long             r_cap;     .. residual capasity
         struct node_st   *head;     .. head node 
         struct arc_st    *sister;   .. opposite arc 
         struct arc_st    *next;     .. next arc with the same tail 
         ..................
       }
   arc;

   type   node   must contain the field 'first': 

   typedef
     struct node_st
       {
          arc_st        *first;    ..  first outgoing arc 
          ....................
       }
    node;
*/

/* ----------------------------------------------------------------- */

int parse
(
/* all parameters are output */
struct Graph *graph,           /* input graph */
long    *n_ad,                 /* address of the number of nodes */
long    *m_ad,                 /* address of the number of arcs */
node    **nodes_ad,            /* address of the array of nodes */
arc     **arcs_ad,             /* address of the array of arcs */
long    **cap_ad,              /* address of the array of capasities */
node    **source_ad,           /* address of the pointer to the source */
node    **sink_ad,             /* address of the pointer to the source */
long    *node_min_ad           /* address of the minimal node */
)
{

#define MAXLINE       100	/* max line length in the input file */
#define ARC_FIELDS      3	/* no of fields in arc line  */
#define NODE_FIELDS     2	/* no of fields in node line  */
#define P_FIELDS        3       /* no of fields in problem line */
#define PROBLEM_TYPE "max"      /* name of problem type*/


long    n,                      /* internal number of nodes */
        node_min,               /* minimal no of node  */
        node_max,               /* maximal no of nodes */
       *arc_first,              /* internal array for holding
                                     - node degree
                                     - position of the first outgoing arc */
       *arc_tail,               /* internal array: tails of the arcs */
        source,                 /* no of the source */
        sink,                   /* no of the sink */
        /* temporary variables carrying no of nodes */
        head, tail, i;

long    m,                      /* internal number of arcs */
        /* temporary variables carrying no of arcs */
        last, arc_num, arc_new_num;

node    *nodes,                 /* pointer to the node structure */
        *head_p,
        *ndp;

arc     *arcs,                  /* pointer to the arc structure */
        *arc_current,
        *arc_new,
        *arc_tmp;

long    *acap,                  /* array of capasities */
        cap;                    /* capasity of the current arc */

long	pos_current=0;          /* 2*no_alines */

int     edge_num;               /* temporary */

/* --------------------------------------------------------------- */

n = graph->node_max - graph->node_min;
m = graph->m;

node_min = 0;
node_max = n;

/* allocating memory for  'nodes', 'arcs'  and internal arrays */
nodes    = (node*) calloc ( n+2, sizeof(node) );
arcs     = (arc*)  calloc ( 2*m+1, sizeof(arc) );
arc_tail = (long*) calloc ( 2*m,   sizeof(long) ); 
arc_first= (long*) calloc ( n+2, sizeof(long) );
acap     = (long*) calloc ( 2*m, sizeof(long) );

/* arc_first [ 0 .. n+1 ] = 0 - initialized by calloc */

if ( nodes == NULL || arcs == NULL || 
 arc_first == NULL || arc_tail == NULL )
/* memory is not allocated */
    { return 1; }
		     
/* setting pointer to the first arc */
arc_current = arcs;

/* reading source and sink */
source = graph->source - graph->node_min;
sink = graph->sink - graph->node_min;

/* --------------------------------------------------------------- */

/* The main loop:
        -  reads the line of the input,
        -  analises its type,
        -  checks correctness of parameters,
        -  puts data to the arrays,
        -  does service functions
*/

graph->current = graph->first;

while ( graph->current )
{
	for ( edge_num = 0; edge_num < graph->current->num; edge_num++ )
	{
		tail = graph->current->edge[edge_num].from - graph->node_min;
		head = graph->current->edge[edge_num].to - graph->node_min;
		cap  = graph->current->edge[edge_num].cap;

		/* no of arcs incident to node i is stored in arc_first[i+1] */
		arc_first[tail + 1] ++; 
		arc_first[head + 1] ++;

                /* storing information about the arc */
		arc_tail[pos_current]        = tail;
		arc_tail[pos_current+1]      = head;
		arc_current       -> head    = nodes + head;
		arc_current       -> r_cap    = cap;
		arc_current       -> sister  = arc_current + 1;
		( arc_current + 1 ) -> head    = nodes + tail;
		( arc_current + 1 ) -> r_cap    = 0;
		( arc_current + 1 ) -> sister  = arc_current;

		arc_current += 2;
		pos_current += 2;
	}
  
	graph->current = graph->current->next;
	free(graph->first);
	graph->first = graph->current;
}     /* end of input loop */

/********** ordering arcs - linear time algorithm ***********/

/* first arc from the first node */
( nodes + node_min ) -> first = arcs;

/* before below loop arc_first[i+1] is the number of arcs outgoing from i;
   after this loop arc_first[i] is the position of the first 
   outgoing from node i arcs after they would be ordered;
   this value is transformed to pointer and written to node.first[i]
   */
 
for ( i = node_min + 1; i <= node_max + 1; i ++ ) 
  {
    arc_first[i]          += arc_first[i-1];
    ( nodes + i ) -> first = arcs + arc_first[i];
  }


for ( i = node_min; i < node_max; i ++ ) /* scanning all the nodes  
                                            exept the last*/
  {

    last = ( ( nodes + i + 1 ) -> first ) - arcs;
                             /* arcs outgoing from i must be cited    
                              from position arc_first[i] to the position
                              equal to initial value of arc_first[i+1]-1  */

    for ( arc_num = arc_first[i]; arc_num < last; arc_num ++ )
      { tail = arc_tail[arc_num];

	while ( tail != i )
          /* the arc no  arc_num  is not in place because arc cited here
             must go out from i;
             we'll put it to its place and continue this process
             until an arc in this position would go out from i */

	  { arc_new_num  = arc_first[tail];
	    arc_current  = arcs + arc_num;
	    arc_new      = arcs + arc_new_num;
	    
	    /* arc_current must be cited in the position arc_new    
	       swapping these arcs:                                 */

	    head_p               = arc_new -> head;
	    arc_new -> head      = arc_current -> head;
	    arc_current -> head  = head_p;

	    cap                 = arc_new -> r_cap;
	    arc_new -> r_cap     = arc_current -> r_cap;
	    arc_current -> r_cap = cap;

	    if ( arc_new != arc_current -> sister )
	      {
	        arc_tmp                = arc_new -> sister;
	        arc_new  -> sister     = arc_current -> sister;
	        arc_current -> sister  = arc_tmp;

                ( arc_current -> sister ) -> sister = arc_current;
		( arc_new     -> sister ) -> sister = arc_new;
	      }

	    arc_tail[arc_num] = arc_tail[arc_new_num];
	    arc_tail[arc_new_num] = tail;

	    /* we increase arc_first[tail]  */
	    arc_first[tail] ++ ;

            tail = arc_tail[arc_num];
	  }
      }
    /* all arcs outgoing from  i  are in place */
  }       

/* -----------------------  arcs are ordered  ------------------------- */

/*----------- constructing lists ---------------*/


  for ( ndp = nodes + node_min; ndp <= nodes + node_max;  ndp ++ )
      ndp -> first = (arc*) NULL;

  for ( arc_current = arcs + (2*m-1); arc_current >= arcs; arc_current -- )
    {
      arc_num = arc_current - arcs;
      tail = arc_tail [arc_num];
      ndp = nodes + tail;
      arc_current -> next = ndp -> first;
      ndp -> first = arc_current;
    }


/* ----------- assigning output values ------------*/
*m_ad = m;
*n_ad = node_max - node_min + 1;
*source_ad = nodes + source;
*sink_ad   = nodes + sink;
*node_min_ad = node_min;
*nodes_ad = nodes + node_min;
*arcs_ad = arcs;
*cap_ad = acap;

for ( arc_current = arcs, arc_num = 0; 
      arc_num < 2*m;
      arc_current ++, arc_num ++
    )
     acap [ arc_num ] = arc_current -> r_cap; 

/* free internal memory */
free ( arc_first ); free ( arc_tail );

/* Thanks God! all is done */
return (0);

}
/* --------------------   end of parser  -------------------*/
