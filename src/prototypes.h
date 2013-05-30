#ifndef DEFINITIONS_H_INCLUDED
#define DEFINITIONS_H_INCLUDED

/* Précision des comparaisons : */

#define EPSILON 1E-10


/* Prototypes des fonctions : */

int comparerDoubles(const void *dble1, const void *dble2);

void calculerPentes(const int *nbPts, const int *loc, const double *x, const double *y,
                     double *pentes, int *ind1, int *ind2, int *nbPentesInfinies);

void indCrsp(const int *nloc, const int *loc1, const int *loc2, const int *crsp,
              int *ind1, int *ind2);

void pondererClasseSp(double ponderation[2], const int coupure,
                     const int nLoc, const double *x, const double *y, const int *n,
                     const int *tabOrd , const int eff1, const int eff2);

void partObl(const int *effMinFils, const int *ponderer,
             const int *nLoc, const double *x, const double *y, const double *z, const int *n,
             double *somCum, int *effCum, const int *nbPentesEgalesMax,
             const double *pentes, const int *ind1, const int *ind2,
             double *regleBin, double *vic, int *partition);

void pondererClasseVql(double ponderation[2],
                       const int *nMod, const int numVql, const int mod,
                       const double *xMod, const double *yMod,
                       const double *x2Mod, const double *y2Mod, const double *xyMod,
                       const int eff1, const int eff2,
                       const int d);

void partVql(const int *nObs, const double *xMod, const double *yMod,
             const double *x2Mod, const double *y2Mod, const double *xyMod,
             const int *nVql, const int *nMod,
             const double *zMod, const double *somzT, const int *effMod,
             const int *ponderer, const int *effMinFils,
             int *numVqt, int *mod, double *vic);

void pondererClasseVqt(double ponderation[2], const int coupure,
                     const int nObs, const int numVqt, const double *x, const double *y,
                     const int *tabOrdre , const int eff1, const int eff2);

void partVqt(const int *nObs, const double *x, const double *y, const double *z, const double *vqt,
             const int *nVqt, const int *ordre, const double *somzT,
             const int *ponderer, const int *effMinFils,
             int *numVqt, double *seuil, double *vic);

void interSegments(const int *m, const double *x, const double *y, const double *coeff,
                   int *coupure, double *coupureX, double *coupureY);

void classeAdj(const int *nc, const int *indC, const double *x, const double *y, int *adj);

void sgmtsGrf(const int *indC, const double *x1, const double *y1, const double *x2, const double *y2,
              double *grf);

#endif /* DEFINITIONS_H_INCLUDED */
