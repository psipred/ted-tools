/* BINARY SORTS ON (string)a WITH POINTERS p (a IS UNCHANGED), BIGGEST TO TOP */

sortstr ( a, p, n, init_pointers )
	char	**a;
	short	*p;
	int	n,
		init_pointers;
{
	register	i, j, m;
	int	pair = 0, floating = 0, integer = 0;

	if ( n <= 0 ) return;							/* DEAL WITH TRIVIAL CASES */
	if ( n == 1 ) {
		p[0] = 0;
		return;
	}
	if( init_pointers ) for ( i = 0; i < n; i++ ) p[i] = i;			/* INITIALISE THE POINTERS */
	m = n;
	while (1) { 								/* LOOP OF DECREASING SWOP SPAN	*/
		m = m/2;
		if(m == 0) return;
		for ( i = 0; i < n-m; i++ ) {					/* LOOP OVER SWOPS */
			j = i;							/* LOOP TO BUBBLE UP */
			while (1) 
			{	register jp;
				if (strcmp(a[p[j+m]],a[p[j]]) >= 0) break;
				jp = p[j];
				p[j] = p[j+m];					/* EXCHANGE POINTERS */
				p[j+m] = jp;
				j = j-m;					/* RETURN TO TOP OF LIST */
				if( j < 0 ) break;				/* ONLY IF ROOM */
			}
		}
	}
}

