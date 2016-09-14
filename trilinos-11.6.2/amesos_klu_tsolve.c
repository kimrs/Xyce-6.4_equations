/* ========================================================================== */
/* === KLU_tsolve =========================================================== */
/* ========================================================================== */

/* Solve A'x=b using the symbolic and numeric objects from KLU_analyze
 * (or KLU_analyze_given) and KLU_factor.  Note that no iterative refinement is
 * performed.  Uses Numeric->Xwork as workspace (undefined on input and output),
 * of size 4n Entry's (note that columns 2 to 4 of Xwork overlap with
 * Numeric->Iwork).
 */

#include "amesos_klu_internal.h"

Int KLU_tsolve
(
    /* inputs, not modified */
    KLU_symbolic *Symbolic,
    KLU_numeric *Numeric,
    Int d,		    /* leading dimension of B */
    Int nrhs,		    /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    double B [ ],	    /* size n*nrhs, in column-oriented form, with
			     * leading dimension d. */
#ifdef COMPLEX
    Int conj_solve,	    /* TRUE for conjugate transpose solve, FALSE for
			     * array transpose solve.  Used for the complex
			     * case only. */
#endif
    /* --------------- */
    KLU_common *Common
)
{
    /* EDIT KIM */
    for(int i = 0; i < 16; ++i) 
        printf("B[%i]: %f\n", i, B[i]);

    printf("tSolve invoked with: \n");
    printf("Symbolic: \n");
    printf("\tsymmetry: %f\n"         ,Symbolic->symmetry);
    printf("\test_flops: %f\n"        ,Symbolic->est_flops);
    printf("\tlnz: %f\n"              ,Symbolic->lnz);
    printf("\tunz: %f\n"              ,Symbolic->unz);
    printf("\tLnz: %d\n"              ,*Symbolic->Lnz);
    printf("\tn: %d\n"                ,Symbolic->n);
    printf("\tnz: %d\n"               ,Symbolic->nz);
    printf("\tP: %d\n"                ,*Symbolic->P);  
    printf("\tQ: %d\n"                ,*Symbolic->Q);
    printf("\tR: %d\n"                ,*Symbolic->R);
    printf("\tnzoff: %d\n"            ,Symbolic->nzoff);
    printf("\tnblocks: %d\n"          ,Symbolic->nblocks);
    printf("\tmaxblock: %d\n"         ,Symbolic->maxblock);
    printf("\tordering: %d\n"         ,Symbolic->ordering);
    printf("\tdo_btf: %d\n"           ,Symbolic->do_btf);
    printf("\tstructural_rank: %d\n"  ,Symbolic->structural_rank);

    printf("Numeric: \n");
    printf("\tn : %d\n", Numeric->n);
    printf("\tnblocks: %d\n", Numeric->nblocks);
    printf("\tlnz: %d\n", Numeric->lnz);
    printf("\tunz: %d\n", Numeric->unz);
    printf("\tmax_lnz_block: %d\n", Numeric->max_lnz_block);
    printf("\tmax_unz_block: %d\n", Numeric->max_unz_block);
    printf("\tPnum: %d\n", *Numeric->Pnum);
    printf("\tPinv: %d\n", *Numeric->Pinv);
    printf("\tLip: %d\n",  *Numeric->Lip);
    printf("\tUip: %d\n",  *Numeric->Uip);
    printf("\tLlen: %d\n", *Numeric->Llen);
    printf("\tUlen: %d\n", *Numeric->Ulen);
    printf("\tUdiag: %d\n", Numeric->Udiag);
    double * Udiag_ = Numeric->Udiag;
    printf("\t\tUdiag[0]: %f\n",  Udiag_[0]);
    printf("\t\tUdiag[1]: %f\n",  Udiag_[1]);
    printf("\t\tUdiag[2]: %f\n",  Udiag_[2]);
    printf("\t\tUdiag[3]: %f\n",  Udiag_[3]);
    printf("\t\tUdiag[4]: %f\n",  Udiag_[4]);
    printf("\t\tUdiag[5]: %f\n",  Udiag_[5]);
    printf("\t\tUdiag[6]: %f\n",  Udiag_[6]);
    printf("\t\tUdiag[7]: %f\n",  Udiag_[7]);
    printf("\t\tUdiag[8]: %f\n",  Udiag_[8]);
    printf("\t\tUdiag[9]: %f\n",  Udiag_[9]);
    printf("\t\tUdiag[10]: %f\n", Udiag_[10]);
    printf("\t\tUdiag[11]: %f\n", Udiag_[11]);
    printf("\tLUsize: %d\n", Numeric->LUsize);
    printf("\tRs: %f\n", Numeric->Rs);
    printf("\tworksize: %d\n", Numeric->worksize);
    printf("\tIwork: %d\n", Numeric->Iwork);
    printf("\tOffp: %d\n", Numeric->Offp);
    printf("\tOffi: %d\n", Numeric->Offi);
    printf("\tOffx: %d\n", Numeric->Offx);
    printf("\tnzoff: %d\n", Numeric->nzoff);

    printf("Common: \n");
    printf("\thalt_if_singular: %d\n", Common->halt_if_singular);
    printf("\tstatus: %d\n", Common->status);
    printf("\tnrealloc: %d\n", Common->nrealloc);
    printf("\tstructural_rank: %d\n", Common->structural_rank);
    printf("\tnumerical_rank: %d\n", Common->numerical_rank);
    printf("\tsingular_col: %d\n", Common->singular_col);
    printf("\tnoffdiag: %d\n", Common->noffdiag);
    printf("\tflops: %f\n", Common->flops);
    printf("\trcond: %f\n", Common->rcond);
    printf("\tcondest: %f\n", Common->condest);
    printf("\trgrowth: %f\n", Common->rgrowth);
    printf("\twork: %f\n", Common->work);
    /* end kim*/

    Entry x [4], offik, s ;
    double rs, *Rs ;
    Entry *Offx, *X, *Bz, *Udiag ;
    Int *Q, *R, *Pnum, *Offp, *Offi, *Lip, *Uip, *Llen, *Ulen ;
    Unit **LUbx ;
    Int k1, k2, nk, k, block, pend, n, p, nblocks, chunk, nr, i ;
#ifdef KLU_ENABLE_OPENMP
    Entry *X1 ;
    int tid ;
#endif

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (Common == NULL)
    {
	return (FALSE) ;
    }
    if (Numeric == NULL || Symbolic == NULL || d < Symbolic->n || nrhs < 0 ||
	B == NULL)
    {
	Common->status = KLU_INVALID ;
	return (FALSE) ;
    }
    Common->status = KLU_OK ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    Bz = (Entry *) B ;
    n = Symbolic->n ;
    nblocks = Symbolic->nblocks ;
    Q = Symbolic->Q ;
    R = Symbolic->R ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Numeric object */
    /* ---------------------------------------------------------------------- */

    ASSERT (nblocks == Numeric->nblocks) ;
    Pnum = Numeric->Pnum ;
    Offp = Numeric->Offp ;
    Offi = Numeric->Offi ;
    Offx = (Entry *) Numeric->Offx ;

    Lip  = Numeric->Lip ;
    Llen = Numeric->Llen ;
    Uip  = Numeric->Uip ;
    Ulen = Numeric->Ulen ;
    LUbx = (Unit **) Numeric->LUbx ;
    Udiag = Numeric->Udiag ;

    Rs = Numeric->Rs ;
#ifdef KLU_ENABLE_OPENMP
    X1 = (Entry *) Numeric->Xwork ;
#else
    X = (Entry *) Numeric->Xwork ;
#endif
    ASSERT (KLU_valid (n, Offp, Offi, Offx)) ;

    /* ---------------------------------------------------------------------- */
    /* solve in chunks of 4 columns at a time */
    /* ---------------------------------------------------------------------- */

#ifdef KLU_ENABLE_OPENMP
#pragma omp parallel for schedule(guided) private(tid, nr, k, block, k1, k2, nk, pend, p, s, i, Bz, X, x, offik, rs) 
#endif
    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
    {
#ifdef KLU_ENABLE_OPENMP
        Bz  = ((Entry *) B) + d*chunk ;
        tid = omp_get_thread_num();
        X = X1 + (tid * n * 4);
#endif


	/* ------------------------------------------------------------------ */
	/* get the size of the current chunk */
	/* ------------------------------------------------------------------ */

	nr = MIN (nrhs - chunk, 4) ;

	/* ------------------------------------------------------------------ */
	/* permute the right hand side, X = Q'*B */
	/* ------------------------------------------------------------------ */

	switch (nr)
	{

	    case 1:

		for (k = 0 ; k < n ; k++)
		{
		    X [k] = Bz  [Q [k]] ;
            printf("X[%d]=Bz(%f)\tQ[%d]= %f\n", k, Bz[Q[k]], k, Q[k]); /*EDIT KIM*/ 
		}
		break ;

	    case 2:

		for (k = 0 ; k < n ; k++)
		{
		    i = Q [k] ;
		    X [2*k    ] = Bz [i      ] ;
		    X [2*k + 1] = Bz [i + d  ] ;
		}
		break ;

	    case 3:

		for (k = 0 ; k < n ; k++)
		{
		    i = Q [k] ;
		    X [3*k    ] = Bz [i      ] ;
		    X [3*k + 1] = Bz [i + d  ] ;
		    X [3*k + 2] = Bz [i + d*2] ;
		}
		break ;

	    case 4:

		for (k = 0 ; k < n ; k++)
		{
		    i = Q [k] ;
		    X [4*k    ] = Bz [i      ] ;
		    X [4*k + 1] = Bz [i + d  ] ;
		    X [4*k + 2] = Bz [i + d*2] ;
		    X [4*k + 3] = Bz [i + d*3] ;
		}
		break ;

	}

	/* ------------------------------------------------------------------ */
	/* solve X = (L*U + Off)'\X */
	/* ------------------------------------------------------------------ */

	for (block = 0 ; block < nblocks ; block++)
	{

	    /* -------------------------------------------------------------- */
	    /* the block of size nk is from rows/columns k1 to k2-1 */
	    /* -------------------------------------------------------------- */

	    k1 = R [block] ;
	    k2 = R [block+1] ;
	    nk = k2 - k1 ;
	    PRINTF (("tsolve %d, k1 %d k2-1 %d nk %d\n", block, k1,k2-1,nk)) ;

	    /* -------------------------------------------------------------- */
	    /* block back-substitution for the off-diagonal-block entries */
	    /* -------------------------------------------------------------- */

	    if (block > 0)
	    {
		switch (nr)
		    {

		    case 1:

			for (k = k1 ; k < k2 ; k++)
			{
			    pend = Offp [k+1] ;
			    for (p = Offp [k] ; p < pend ; p++)
			    {
#ifdef COMPLEX
				if (conj_solve)
				{
				    MULT_SUB_CONJ (X [k], X [Offi [p]],
					    Offx [p]) ;
				}
				else
#endif
				{
				    MULT_SUB (X [k], Offx [p], X [Offi [p]]) ;
				}
			    }
			}
			break ;

		    case 2:

			for (k = k1 ; k < k2 ; k++)
			{
			    pend = Offp [k+1] ;
			    x [0] = X [2*k    ] ;
			    x [1] = X [2*k + 1] ;
			    for (p = Offp [k] ; p < pend ; p++)
			    {
				i = Offi [p] ;
#ifdef COMPLEX
				if (conj_solve)
				{
				    CONJ (offik, Offx [p]) ;
				}
				else
#endif
				{
				    offik = Offx [p] ;
				}
				MULT_SUB (x [0], offik, X [2*i]) ;
				MULT_SUB (x [1], offik, X [2*i + 1]) ;
			    }
			    X [2*k    ] = x [0] ;
			    X [2*k + 1] = x [1] ;
			}
			break ;

		    case 3:

			for (k = k1 ; k < k2 ; k++)
			{
			    pend = Offp [k+1] ;
			    x [0] = X [3*k    ] ;
			    x [1] = X [3*k + 1] ;
			    x [2] = X [3*k + 2] ;
			    for (p = Offp [k] ; p < pend ; p++)
			    {
				i = Offi [p] ;
#ifdef COMPLEX
				if (conj_solve)
				{
				    CONJ (offik, Offx [p]) ;
				}
				else
#endif
				{
				    offik = Offx [p] ;
				}
				MULT_SUB (x [0], offik, X [3*i]) ;
				MULT_SUB (x [1], offik, X [3*i + 1]) ;
				MULT_SUB (x [2], offik, X [3*i + 2]) ;
			    }
			    X [3*k    ] = x [0] ;
			    X [3*k + 1] = x [1] ;
			    X [3*k + 2] = x [2] ;
			}
			break ;

		    case 4:

			for (k = k1 ; k < k2 ; k++)
			{
			    pend = Offp [k+1] ;
			    x [0] = X [4*k    ] ;
			    x [1] = X [4*k + 1] ;
			    x [2] = X [4*k + 2] ;
			    x [3] = X [4*k + 3] ;
			    for (p = Offp [k] ; p < pend ; p++)
			    {
				i = Offi [p] ;
#ifdef COMPLEX
				if (conj_solve)
				{
				    CONJ(offik, Offx [p]) ;
				}
				else
#endif
				{
				    offik = Offx [p] ;
				}
				MULT_SUB (x [0], offik, X [4*i]) ;
				MULT_SUB (x [1], offik, X [4*i + 1]) ;
				MULT_SUB (x [2], offik, X [4*i + 2]) ;
				MULT_SUB (x [3], offik, X [4*i + 3]) ;
			    }
			    X [4*k    ] = x [0] ;
			    X [4*k + 1] = x [1] ;
			    X [4*k + 2] = x [2] ;
			    X [4*k + 3] = x [3] ;
			}
			break ;
		    }
	    }

	    /* -------------------------------------------------------------- */
	    /* solve the block system */
	    /* -------------------------------------------------------------- */

	    if (nk == 1)
	    {
#ifdef COMPLEX
		if (conj_solve)
		{
		    CONJ (s, Udiag [k1]) ;
		}
		else
#endif
		{
		    s = Udiag [k1] ;
		}
		switch (nr)
		{

		    case 1:
			DIV (X [k1], X [k1], s) ;
			break ;

		    case 2:
			DIV (X [2*k1], X [2*k1], s) ;
			DIV (X [2*k1 + 1], X [2*k1 + 1], s) ;
			break ;

		    case 3:
			DIV (X [3*k1], X [3*k1], s) ;
			DIV (X [3*k1 + 1], X [3*k1 + 1], s) ;
			DIV (X [3*k1 + 2], X [3*k1 + 2], s) ;
			break ;

		    case 4:
			DIV (X [4*k1], X [4*k1], s) ;
			DIV (X [4*k1 + 1], X [4*k1 + 1], s) ;
			DIV (X [4*k1 + 2], X [4*k1 + 2], s) ;
			DIV (X [4*k1 + 3], X [4*k1 + 3], s) ;
			break ;

		}
	    }
	    else
	    {
    /* EDIT KIM */
    printf("before solve:\n");
    for(int i = 0; i < 16; ++i) 
        printf("Udiag[%d]: %f\n", i, Udiag[i]); 

		KLU_utsolve (nk, Uip + k1, Ulen + k1, LUbx [block],
			Udiag + k1, nr, /* note Udiag is slightly diferent */
#ifdef COMPLEX
			conj_solve,
#endif
			X + nr*k1) ;
        /* EDIT KIM */
        printf("out of utsolve:\n");
        for(int i = 0; i < 16; ++i) 
            printf("X[%d]: %e\n", i, X[i]); 

		KLU_ltsolve (nk, Lip + k1, Llen + k1, LUbx [block], nr, /*NOTE Lip is different here than from amix*/
#ifdef COMPLEX
			conj_solve,
#endif
			X + nr*k1) ;
         /* EDIT KIM */
         printf("out of ltsolve:\n");
         for(int i = 0; i < 16; ++i) 
             printf("X[%d]: %e\n", i, X[i]); 
	    }
	}

	/* ------------------------------------------------------------------ */
	/* scale and permute the result, Bz  = P'(R\X) */
	/* ------------------------------------------------------------------ */

	if (Rs == NULL)
	{

	    /* no scaling */
	    switch (nr)
	    {

		case 1:

		    for (k = 0 ; k < n ; k++) {
			    Bz  [Pnum [k]] = X [k] ;
                /*EDIT KIM */
                printf("Bz[Pnum[%d]] = X[%d]\n", k, k);
                printf("Bz[%d] = %g\n", Pnum[k], X[k]);
		    }
            /* EDIT KIM */
            for(int i = 0; i < 16; ++i) 
                printf("B[%i]: %e\n", i, B[i]);
		    break ;

		case 2:

		    for (k = 0 ; k < n ; k++)
		    {
			i = Pnum [k] ;
			Bz  [i      ] = X [2*k    ] ;
			Bz  [i + d  ] = X [2*k + 1] ;
		    }
		    break ;

		case 3:

		    for (k = 0 ; k < n ; k++)
		    {
			i = Pnum [k] ;
			Bz  [i      ] = X [3*k    ] ;
			Bz  [i + d  ] = X [3*k + 1] ;
			Bz  [i + d*2] = X [3*k + 2] ;
		    }
		    break ;

		case 4:

		    for (k = 0 ; k < n ; k++)
		    {
			i = Pnum [k] ;
			Bz  [i      ] = X [4*k    ] ;
			Bz  [i + d  ] = X [4*k + 1] ;
			Bz  [i + d*2] = X [4*k + 2] ;
			Bz  [i + d*3] = X [4*k + 3] ;
		    }
		    break ;
	    }

	}
	else
	{

	    switch (nr)
	    {

		case 1:

		    for (k = 0 ; k < n ; k++)
		    {
			SCALE_DIV_ASSIGN (Bz [Pnum [k]], X [k], Rs [k]) ;
		    }
		    break ;

		case 2:

		    for (k = 0 ; k < n ; k++)
		    {
			i = Pnum [k] ;
			rs = Rs [k] ;
			SCALE_DIV_ASSIGN (Bz [i], X [2*k], rs) ;
			SCALE_DIV_ASSIGN (Bz [i + d], X [2*k + 1], rs) ;
		    }
		    break ;

		case 3:

		    for (k = 0 ; k < n ; k++)
		    {
			i = Pnum [k] ;
			rs = Rs [k] ;
			SCALE_DIV_ASSIGN (Bz [i], X [3*k], rs) ;
			SCALE_DIV_ASSIGN (Bz [i + d], X [3*k + 1], rs) ;
			SCALE_DIV_ASSIGN (Bz [i + d*2], X [3*k + 2], rs) ;
		    }
		    break ;

		case 4:

		    for (k = 0 ; k < n ; k++)
		    {
			i = Pnum [k] ;
			rs = Rs [k] ;
			SCALE_DIV_ASSIGN (Bz [i], X [4*k], rs) ;
			SCALE_DIV_ASSIGN (Bz [i + d], X [4*k + 1], rs) ;
			SCALE_DIV_ASSIGN (Bz [i + d*2], X [4*k + 2], rs) ;
			SCALE_DIV_ASSIGN (Bz [i + d*3], X [4*k + 3], rs) ;
		    }
		    break ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* go to the next chunk of B */
	/* ------------------------------------------------------------------ */

#ifndef KLU_ENABLE_OPENMP
	Bz  += d*4 ;
#endif
    }
    /*EDIT KIM */
    for(int i = 0; i < 16; ++i) {
        printf("X[%d]: %f\n", X[i]);
    }
    return (TRUE) ;
}
