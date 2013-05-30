#include <R.h>
#include <math.h>
#include <stdlib.h>
#include "prototypes.h"

void partObl(const int *effMinFils, const int *ponderer,
             const int *nLoc, const double *x, const double *y, const double *z, const int *n,
             double *somCum, int *effCum, const int *nbPentesEgalesMax,
             const double *pentes, const int *ind1, const int *ind2,
             double *coeff, double *vic, int *partition)
{
    /* Variables locales : */
    int i, booleen;
    int nbPerm, p;
    int i1, i2, pente;
    int nbIndActu, indTemp, pos, pos1, pos2;
    double vicTemp, comp, pond[2], pondMax, xTemp, yTemp, yLim;

    /* Tableaux dynamiques : */
    int *ordre, *position, *indActu;


    /* INITIALISATIONS : */

    booleen = 0;
    nbPerm = *nLoc*(*nLoc-1)/2;
    i1 = -1;
    i2 = -1;
    pente = -1;
    nbIndActu = 0;
    indTemp = -1;
    pos = -1;
    pos1 = -1;
    pos2 = -1;
    vicTemp = 0;
    *vic = 0;
    comp = 0;
    pond[0] = 1;
    pond[1] = 1;
    pondMax = 0;
    xTemp = 0;
    yTemp = 0;
    yLim = 0;

    ordre = malloc(*nLoc*sizeof(int));
    position = malloc(*nLoc*sizeof(int));
    indActu = malloc(*nbPentesEgalesMax*sizeof(int));
    if (ordre != NULL  &&  position != NULL  &&  indActu != NULL)
    {
            for (i=0; i<*nLoc; i++)
            {
                ordre[i]=i;
                position[i]=i;
            }
            for (i=0; i<*nbPentesEgalesMax; i++)
                indActu[i] = -1;


            /* RECHERCHE DE LA MEILLEURE PARTITION : */

            for (p=0; p<nbPerm; p++)
            {
                pos1 = position[ind1[p]];
                pos2 = position[ind2[p]];
                indTemp = ordre[pos1];
                ordre[pos1] = ordre[pos2];
                ordre[pos2] = indTemp;
                if (pos1 < pos2)
                {
                    effCum[pos1] += n[ind2[p]] - n[ind1[p]];
                    somCum[pos1] += z[ind2[p]] - z[ind1[p]];
                    pos = pos1;
                }
                else
                {
                    effCum[pos2] += n[ind1[p]] - n[ind2[p]];
                    somCum[pos2] += z[ind1[p]] - z[ind2[p]];
                    pos = pos2;
                }
                indTemp = position[ind1[p]];
                position[ind1[p]] = position[ind2[p]];
                position[ind2[p]] = indTemp;

                booleen = 0;
                for (i=0; i<=nbIndActu; i++)
                    if (pos == indActu[i])
                        booleen = 1;
                if (booleen == 0)
                {
                    indActu[nbIndActu] = pos;
                    nbIndActu++;
                }

                if ((p != nbPerm-1  &&  comparerDoubles(&pentes[p], &pentes[p+1]) != 0)  ||  p == nbPerm-1)
                {
                    for ( i=0; i<nbIndActu; i++)
                    {
                        if (effCum[indActu[i]] < *effMinFils  ||  effCum[*nLoc-1] - effCum[indActu[i]] < *effMinFils)
                            continue;

                        if (*ponderer == 1)
                            pondererClasseSp(pond, indActu[i], *nLoc, x, y, n, ordre, effCum[indActu[i]], effCum[*nLoc-1]-effCum[indActu[i]]);

                        vicTemp = pond[0] *            effCum[indActu[i]]        * (           somCum[indActu[i]]       /         effCum[indActu[i]]          - somCum[*nLoc-1]/effCum[*nLoc-1])
                                                                                 * (           somCum[indActu[i]]       /         effCum[indActu[i]]          - somCum[*nLoc-1]/effCum[*nLoc-1])
                                + pond[1] * (effCum[*nLoc-1]-effCum[indActu[i]]) * ((somCum[*nLoc-1]-somCum[indActu[i]])/(effCum[*nLoc-1]-effCum[indActu[i]]) - somCum[*nLoc-1]/effCum[*nLoc-1])
                                                                                 * ((somCum[*nLoc-1]-somCum[indActu[i]])/(effCum[*nLoc-1]-effCum[indActu[i]]) - somCum[*nLoc-1]/effCum[*nLoc-1]);
                        if (comparerDoubles(&vicTemp, vic) == 1)
                        {
                            if (*ponderer == 1)
                            {
                                if (comparerDoubles(&pond[0], &pond[1]) == 1)
                                    pondMax = pond[0];
                                else
                                    pondMax = pond[1];
                            }
                            i1 = ordre[indActu[i]];
                            i2 = ordre[indActu[i]+1];
                            pente = p;
                            *vic = vicTemp;
                        }
                        else if (*ponderer == 1  &&  i1 != -1  &&  comparerDoubles(&vicTemp, vic) == 0)
                        {
                            if (comparerDoubles(&pond[0], &pond[1]) == 1)
                                comp = pond[0];
                            else
                                comp = pond[1];
                            if (comparerDoubles(&comp, &pondMax) == 1)
                            {
                                pondMax = comp;
                                i1 = ordre[indActu[i]];
                                i2 = ordre[indActu[i]+1];
                                pente = p;
                            }
                        }
                    }
                    nbIndActu = 0;
                }
            }
            free(indActu);
            free(ordre);
            free(position);


            /* STOCKAGE DU RESULTAT : */

            if (pente == -1)
                return;

            xTemp = (x[i1] + x[i2]) / 2;
            yTemp = (y[i1] + y[i2]) / 2;

            if (pente != nbPerm - 1)
            {
                if (pentes[pente+1] > pentes[0])
                    coeff[0] = (pentes[pente] + pentes[pente+1]) / 2;
                else
                    coeff[0] = pentes[pente] * 2;
            }
            else
            {
                if (pentes[pente] > pentes[0])
                    coeff[0] = pentes[pente] * 2;
                else
                    coeff[0] = pentes[pente];
            }

            coeff[1] = yTemp - coeff[0] * xTemp;

            for (i=0; i<*nLoc; i++)
            {
                yLim = coeff[0] * x[i] + coeff[1];
                if (comparerDoubles(&y[i], &yLim) != 1)
                    partition[i] = -1;
                else
                    partition[i] = 1;
            }

            return;
    }

}

void pondererClasseSp(double ponderation[2], const int coupure,
                     const int nbLoc, const double *x, const double *y, const int *n,
                     const int *tabOrdre , const int eff1, const int eff2)
{
    int i;
    double moyX, moyY, varX, varY, covXY, delta;

    moyX = 0;
    moyY = 0;
    varX = 0;
    varY = 0;
    covXY = 0;
    for (i=0; i<=coupure; i++)
    {
        moyX += x[tabOrdre[i]] * n[tabOrdre[i]];
        moyY += y[tabOrdre[i]] * n[tabOrdre[i]];
        varX += x[tabOrdre[i]] * x[tabOrdre[i]] * n[tabOrdre[i]];
        varY += y[tabOrdre[i]] * y[tabOrdre[i]] * n[tabOrdre[i]];
        covXY += x[tabOrdre[i]] * y[tabOrdre[i]] * n[tabOrdre[i]];
    }
    moyX /= eff1;
    moyY /= eff1;
    varX = varX / eff1 - moyX * moyX;
    varY = varY / eff1 - moyY * moyY;
    covXY = covXY / eff1 - moyX * moyY;
    delta = varX * varY - covXY * covXY + eff1;
    ponderation[0] = exp(eff1 / delta) / (1 + exp(eff1 / delta));

    moyX = 0;
    moyY = 0;
    varX = 0;
    varY = 0;
    covXY = 0;
    for (i=coupure+1; i<nbLoc; i++)
    {
        moyX += x[tabOrdre[i]] * n[tabOrdre[i]];
        moyY += y[tabOrdre[i]] * n[tabOrdre[i]];
        varX += x[tabOrdre[i]] * x[tabOrdre[i]] * n[tabOrdre[i]];
        varY += y[tabOrdre[i]] * y[tabOrdre[i]] * n[tabOrdre[i]];
        covXY += x[tabOrdre[i]] * y[tabOrdre[i]] * n[tabOrdre[i]];
    }
    moyX /= eff2;
    moyY /= eff2;
    varX = varX / eff2 - moyX * moyX;
    varY = varY / eff2 - moyY * moyY;
    covXY = covXY / eff2 - moyX * moyY;
    delta = varX * varY - covXY * covXY + eff2;
    ponderation[1] = exp(eff2 / delta) / (1 + exp(eff2 / delta));

    return;
}

void calculerPentes(const int *nbPts, const int *loc, const double *x, const double *y,
                     double *pentes, int *ind1, int *ind2, int *nbPentesInfinies)
{
    int i, j, booleen, compte, nbPentesAChanger, *pentesAChanger;
    double penteMax;

    *nbPentesInfinies = 0;
    penteMax = 0;
    booleen = 0;
    for (i=0; i<*nbPts-1; i++)
        for(j=i+1; j<*nbPts; j++)
        {
            if (comparerDoubles(&x[i], &x[j]) == 0)
                *nbPentesInfinies += 1;
            else
                if (booleen == 0)
                {
                    penteMax = (y[i] - y[j]) / (x[i] - x[j]);
                    booleen = 1;
                }
        }


    pentesAChanger = malloc(*nbPentesInfinies*sizeof(int));
    if (pentesAChanger != NULL)
    {
        compte = 0;
        nbPentesAChanger=0;
        for (i=0; i<*nbPts-1; i++)
            for (j=i+1; j<*nbPts; j++)
            {
                if (comparerDoubles(&x[i], &x[j]) != 0)
                {
                    pentes[compte] = (y[i] - y[j]) / (x[i] - x[j]);
                    if (comparerDoubles(&pentes[compte], &penteMax) == 1)
                        penteMax = pentes[compte];
                }
                else
                {
                    pentesAChanger[nbPentesAChanger] = compte;
                    nbPentesAChanger++;
                }
                ind1[compte] = loc[i];
                ind2[compte] = loc[j];
                compte++;
            }

        for (i=0; i<nbPentesAChanger; i++)
            pentes[pentesAChanger[i]] = penteMax + 1;


        free(pentesAChanger);

        return;
    }

}

/*
    nc : nb de classes
    indC[c] : premier point (indice de x et y) appartenant à la classe c
    x : abscisse d'un point
    y : ordonnee d'un point

    adj : matrice d'adjacence (resultat)
*/

void classeAdj(const int *nc, const int *indC, const double *x, const double *y, int *adj)
{
    int c1, c2, i1, i2, j1, j2, booleen;
    double d1, d2, mi, Mi, mj, Mj;

    for (c1=0; c1<*nc-1; c1++)
        for (c2=c1+1; c2<*nc; c2++)
        {
            booleen = 0;
            for (i1=indC[c1]; i1<indC[c1+1]; i1++)
            {
                if (booleen == 1)
                    break;

                for (j1=indC[c2]; j1<indC[c2+1]; j1++)
                {
                    if (i1 != indC[c1+1]-1)
                        i2 = i1 + 1;
                    else
                        i2 = indC[c1];
                    if (j1 != indC[c2+1]-1)
                        j2 = j1 + 1;
                    else
                        j2 = indC[c2];

                    if (comparerDoubles(&x[i1], &x[i2]) == 0  ||  comparerDoubles(&x[j1], &x[j2]) == 0)
                        continue; /* Un des deux segments est un segment vertical */

                    d1 = (x[i1] - x[i2]) * (y[j1] - y[j2]);
                    d2 = (x[j1] - x[j2]) * (y[i1] - y[i2]);

                    if (comparerDoubles(&d1, &d2) != 0)
                        continue; /* Les 2 segments n'ont pas la meme pente */

                    d1 = (x[i1] - x[i2]) * (x[j1]*y[j2] - x[j2]*y[j1]);
                    d2 = (x[j1] - x[j2]) * (x[i1]*y[i2] - x[i2]*y[i1]);

                    if (comparerDoubles(&d1, &d2) != 0)
                        continue; /* Les 2 droites portant les segments n'ont pas la meme ordonnee a l'origine */

                    /* Maintenant il est certain que  les 2 segments sont portes par la même droite */

                    if (comparerDoubles(&x[i1], &x[i2]) == -1)
                    {
                        mi = x[i1];
                        Mi = x[i2];
                    }
                    else
                    {
                        mi = x[i2];
                        Mi = x[i1];
                    }

                    if (comparerDoubles(&x[j1], &x[j2]) == -1)
                    {
                        mj = x[j1];
                        Mj = x[j2];
                    }
                    else
                    {
                        mj = x[j2];
                        Mj = x[j1];
                    }

                    if (comparerDoubles(&Mi, &mj) == -1  ||  comparerDoubles(&mi, &Mj) == 1)
                        continue; /* Les 2 segments ne ce chevauchent pas */

                    adj[*nc*c1+c2] = 1;
                    adj[*nc*c2+c1] = 1;

                    booleen = 1;
                    break;
                }
            }
        }
    return;
}

int comparerDoubles(const void *dble1, const void *dble2)
{
    int comp;
    double a, b;
    a = 0;
    b = 0;

    if (*(double*)dble1 < *(double*)dble2)
    {
        a = *(double*)dble1;
        b = *(double*)dble2;
        comp = -1;
    }
    else if (*(double*)dble1 > *(double*)dble2)
    {
        a = *(double*)dble2;
        b = *(double*)dble1;
        comp = 1;
    }
    else
        comp = 0;


    if (comp != 0)
    {
        if (b - a < EPSILON)
            comp = 0;
        else
            if ((b - a) < EPSILON*fabs(a))
                comp = 0;
    }
    return comp;
}

void partVql(const int *nObs, const double *xMod, const double *yMod,
             const double *x2Mod, const double *y2Mod, const double *xyMod,
             const int *nVql, const int *nMod,
             const double *zMod, const double *somzT, const int *effMod,
             const int *ponderer, const int *effMinFils,
             int *numVql, int *mod, double *vic)
{
    int i, j, k, d, iTemp, jTemp, effCum;
    double somCum, pond[2], pondMax, comp, vicTemp;

    iTemp = -1;
    jTemp = -1;
    pond[0] = 1;
    pond[1] = 1;
    pondMax = 0;
    comp = 0;
    vicTemp = 0;

    for (i=0; i<*nVql; i++)
    {
        if (i == 0)
            d = 0;
        else
            d = nMod[i-1];

        for (j=d; j<d+nMod[i]; j++)
        {
            effCum = *nObs - effMod[j];
            if (effMod[j] < *effMinFils  ||  effCum < *effMinFils)
                continue;

            somCum = 0;

            for (k=d; k<d+nMod[i]; k++)
                if (k != j)
                {
                    somCum += zMod[k];
                }

            if (*ponderer == 1)
                pondererClasseVql(pond, nMod, i, j, xMod, yMod, x2Mod, y2Mod, xyMod, effMod[j], effCum, d);

            vicTemp = pond[0] * effMod[j] * (zMod[j] / effMod[j] - *somzT / *nObs) * (zMod[j] / effMod[j] - *somzT / *nObs)
                    + pond[1] *  effCum   * (somCum  /  effCum   - *somzT / *nObs) * (somCum  /  effCum   - *somzT / *nObs);

            if (comparerDoubles(vic, &vicTemp) == -1)
            {
                if (*ponderer == 1)
                {
                    if (comparerDoubles(&pond[0], &pond[1]) == 1)
                        pondMax = pond[0];
                    else
                        pondMax = pond[1];
                }
                *vic = vicTemp;
                iTemp = i;
                jTemp = j - d;
            }
            else if (*ponderer == 1  &&  iTemp != -1  &&  comparerDoubles(vic, &vicTemp) == 0)
            {
                if (comparerDoubles(&pond[0], &pond[1]) == 1)
                    comp = pond[0];
                else
                    comp = pond[1];

                if (comparerDoubles(&comp, &pondMax) == 1)
                {
                    pondMax = comp;
                    iTemp = i;
                    jTemp = j - d;
                }
            }
        }
    }
    if (iTemp == -1)
        return;

    *numVql = iTemp + 1;
    *mod = jTemp + 1;

    return;
}

void pondererClasseVql(double ponderation[2],
                       const int *nMod, const int numVql, const int mod,
                       const double *xMod, const double *yMod,
                       const double *x2Mod, const double *y2Mod, const double *xyMod,
                       const int eff1, const int eff2,
                       const int d)
{
    int i;
    double moyX, moyY, varX, varY, covXY, delta;

    moyX = 0;
    moyY = 0;
    varX = 0;
    varY = 0;
    covXY = 0;
    delta = 0;


    moyX = xMod[mod] / eff1;
    moyY = yMod[mod] / eff1;
    varX = x2Mod[mod] / eff1 - moyX;
    varY = y2Mod[mod] / eff1 - moyY;
    covXY = xyMod[mod] / eff1 - moyX * moyY;
    delta = varX * varY - covXY * covXY + eff1;
    ponderation[0] = exp(eff1 / delta) / (1 + exp(eff1 / delta));

    moyX = 0;
    moyY = 0;
    varX = 0;
    varY = 0;
    covXY = 0;
    delta = 0;

    for (i=d; i<d+nMod[numVql]; i++)
    {
        if (i != mod)
        {
            moyX += xMod[i];
            moyY += yMod[i];
            varX += x2Mod[i];
            varY += y2Mod[i];
            covXY += xyMod[i];
        }
    }
    moyX /= eff2;
    moyY /= eff2;
    varX = varX / eff2 - moyX * moyX;
    varY = varY / eff2 - moyY * moyY;
    covXY = covXY / eff2 - moyX * moyY;
    delta = varX * varY - covXY * covXY + eff2;
    ponderation[1] = exp(eff2 / delta) / (1 + exp(eff2 / delta));

    return;
}

void partVqt(const int *nObs, const double *x, const double *y, const double *z, const double *vqt,
             const int *nVqt, const int *ordre, const double *somzT,
             const int *ponderer, const int *effMinFils,
             int *numVqt, double *seuil, double *vic)
{
    /* Variables : */
    int i, j, iTemp, jTemp, effCum;
    double somCum, vicTemp, comp, pond[2], pondMax;


    /* Balayage des decoupages : */

    vicTemp = 0;
    somCum = 0;
    effCum = 0;
    comp = 0;
    pond[0] = 1;
    pond[1] = 1;
    pondMax = 0;
    iTemp = -1;
    jTemp = -1;
    for (i=0; i<*nVqt; i++)
    {
        effCum = 0;
        somCum = 0;
        for (j=0; j<*nObs-1; j++)
        {
            effCum ++;
            somCum += z[ordre[*nObs * i + j]];

            if (effCum < *effMinFils  ||  *nObs - effCum < *effMinFils)
                continue;

            if (comparerDoubles(&vqt[*nObs * i + ordre[*nObs * i + j]], &vqt[*nObs * i + ordre[*nObs * i + j + 1]]) != 0)
            {
                if (*ponderer == 1)
                    pondererClasseVqt(pond, *nObs * i + j, *nObs, i, x, y, ordre, effCum, *nObs-effCum);

                vicTemp = pond[0] *      effCum      * (      somCum      /      effCum      - *somzT / *nObs)
                                                     * (      somCum      /      effCum      - *somzT / *nObs)
                        + pond[1] * (*nObs - effCum) * ((*somzT - somCum) / (*nObs - effCum) - *somzT / *nObs)
                                                     * ((*somzT - somCum) / (*nObs - effCum) - *somzT / *nObs);

                if (comparerDoubles(vic, &vicTemp) == -1)
                {
                    if (*ponderer == 1)
                    {
                        if (comparerDoubles(&pond[0], &pond[1]) == 1)
                            pondMax = pond[0];
                        else
                            pondMax = pond[1];
                    }
                    *vic = vicTemp;
                    iTemp = i;
                    jTemp = j;
                }
                else if (*ponderer == 1  &&  iTemp != -1  &&  comparerDoubles(vic, &vicTemp) == 0)
                {
                    if (comparerDoubles(&pond[0], &pond[1]) == 1)
                        comp = pond[0];
                    else
                        comp = pond[1];

                    if (comparerDoubles(&comp, &pondMax) == 1)
                    {
                        pondMax = comp;
                        iTemp = i;
                        jTemp = j;
                    }
                }
            }
        }
    }
    if (iTemp == -1)
        return;

    *numVqt = iTemp + 1;
    *seuil = (vqt[*nObs * iTemp + ordre[*nObs * iTemp + jTemp]] + vqt[*nObs * iTemp + ordre[*nObs * iTemp + jTemp]]) / 2;

    return;
}

void pondererClasseVqt(double ponderation[2], const int coupure,
                     const int nObs, const int numVqt, const double *x, const double *y,
                     const int *tabOrdre , const int eff1, const int eff2)
{
    int i;
    double moyX, moyY, varX, varY, covXY, delta;

    moyX = 0;
    moyY = 0;
    varX = 0;
    varY = 0;
    covXY = 0;
    for (i=nObs*numVqt; i<=coupure; i++)
    {
        moyX += x[tabOrdre[i]];
        moyY += y[tabOrdre[i]];
        varX += x[tabOrdre[i]] * x[tabOrdre[i]];
        varY += y[tabOrdre[i]] * y[tabOrdre[i]];
        covXY += x[tabOrdre[i]] * y[tabOrdre[i]];
    }
    moyX /= eff1;
    moyY /= eff1;
    varX = varX / eff1 - moyX * moyX;
    varY = varY / eff1 - moyY * moyY;
    covXY = covXY / eff1 - moyX * moyY;
    delta = varX * varY - covXY * covXY + eff1;
    ponderation[0] = exp(eff1 / delta) / (1 + exp(eff1 / delta));

    moyX = 0;
    moyY = 0;
    varX = 0;
    varY = 0;
    covXY = 0;
    for (i=coupure+1; i<nObs*(numVqt+1); i++)
    {
        moyX += x[tabOrdre[i]];
        moyY += y[tabOrdre[i]];
        varX += x[tabOrdre[i]] * x[tabOrdre[i]];
        varY += y[tabOrdre[i]] * y[tabOrdre[i]];
        covXY += x[tabOrdre[i]] * y[tabOrdre[i]];
    }
    moyX /= eff2;
    moyY /= eff2;
    varX = varX / eff2 - moyX * moyX;
    varY = varY / eff2 - moyY * moyY;
    covXY = covXY / eff2 - moyX * moyY;
    delta = varX * varY - covXY * covXY + eff2;
    ponderation[1] = exp(eff2 / delta) / (1 + exp(eff2 / delta));

    return;
}

void indCrsp(const int *nLoc, const int *crsp, const int *loc1, const int *loc2,
             int *ind1, int *ind2)
{
    int i, j, nPerm, bool1, bool2;

    nPerm = *nLoc*(*nLoc-1)/2;

    for (i=0; i<nPerm; i++)
    {
        bool1 = 0;
        bool2 = 0;
        for (j=0; j<*nLoc; j++)
        {
            if (loc1[i] == crsp[j])
            {
                ind1[i] = j;
                bool1 = 1;
            }
            if (loc2[i] == crsp[*nLoc-1-j])
            {
                ind2[i] = *nLoc-1-j;
                bool2 = 1;
            }
            if (bool1 == 1  &&  bool2 == 1)
                break;
        }
    }
    return;
}

void interSegments(const int *m, const double *x, const double *y, const double *coeff,
                   int *coupure, double *coupureX, double *coupureY)
{
    int i, j, compte;
    double minX, maxX, minY, maxY, a, b, u, v;

    for (i=0; i<*m; i++)
        coupure[i] = 3;


    compte = 0;
    for (i=0; i<*m; i++)
    {
        if (compte == 0)
            coupure[i] = 1;
        else if (compte == 1)
            coupure[i] = 2;
        else
            break;

        if (i != *m-1)
            j=i+1;
        else
            j=0;

        if (comparerDoubles(&x[i], &x[j]) == 0  &&  comparerDoubles(&y[i], &y[j]) == 0)
            continue;

        if (comparerDoubles(&y[i], &y[j]) == 1)
        {
            minY = y[j];
            maxY = y[i];
        }
        else
        {
            minY = y[i];
            maxY = y[j];
        }

        if (comparerDoubles(&x[i], &x[j]) == 0)
        {
            u = x[i];
            v = coeff[0] * u + coeff[1];

            if (comparerDoubles(&v, &minY) != -1  &&  comparerDoubles(&v, &maxY) != 1)
            {
                if (compte == 1  &&  comparerDoubles(&u, &coupureX[0]) == 0  &&  comparerDoubles(&v, &coupureY[0]) == 0)
                        continue;
                coupureX[compte] = u;
                coupureY[compte] = v;
                compte++;
            }
        }
        else
        {
            if (comparerDoubles(&x[i], &x[j]) == 1)
            {
                minX = x[j];
                maxX = x[i];
            }
            else
            {
                minX = x[i];
                maxX = x[j];
            }

            a = (y[i] - y[j]) / (x[i] - x[j]);
            b = y[i] - a * x[i];

            if (comparerDoubles(&a, &coeff[0]) != 0)
            {
                u = (b - coeff[1]) / (coeff[0] - a);
                v = coeff[0] * u + coeff[1];
            }
            else
                continue;

            if (comparerDoubles(&u, &minX) != -1  &&  comparerDoubles(&u, &maxX) != 1
            &&  comparerDoubles(&v, &minY) != -1  &&  comparerDoubles(&v, &maxY) != 1)
            {
                if (compte == 1  &&  comparerDoubles(&u, &coupureX[0]) == 0  &&  comparerDoubles(&v, &coupureY[0]) == 0)
                {
                    continue;
                }
                coupureX[compte] = u;
                coupureY[compte] = v;
                compte++;
            }
        }
    }
}

/*
    indC[c] : dernier point (indice de x et y) appartenant à la classe c
    x : abscisse d'un point
    y : ordonnee d'un point

    grf : segments communs à la classe 1 et la classe 2 (résultat)
*/

void sgmtsGrf(const int *indC, const double *x1, const double *y1, const double *x2, const double *y2,
              double *grf)
{
    int i, j, compte;
    double d1, d2, mxi, myi, Mxi, Myi, mxj, myj, Mxj, Myj;

    compte = 0;

    for (i=0; i<indC[0]; i++)
        for (j=indC[0]; j<indC[1]; j++)
        {
            if (comparerDoubles(&x1[i], &x2[i]) == 0  ||  comparerDoubles(&x1[j], &x2[j]) == 0)
                continue; /* Un des deux segments est un segment vertical */

            d1 = (x1[i] - x2[i]) * (y1[j] - y2[j]);
            d2 = (x1[j] - x2[j]) * (y1[i] - y2[i]);

            if (comparerDoubles(&d1, &d2) != 0)
                continue; /* Les 2 segments n'ont pas la meme pente */

            d1 = (x1[i] - x2[i]) * (x1[j]*y2[j] - x2[j]*y1[j]);
            d2 = (x1[j] - x2[j]) * (x1[i]*y2[i] - x2[i]*y1[i]);

            if (comparerDoubles(&d1, &d2) != 0)
                continue; /* Les 2 droites portant les segments n'ont pas la meme ordonnee a l'origine */


            /* Maintenant il est certain que  les 2 segments sont portes par la même droite */

            if (comparerDoubles(&x1[i], &x2[i]) == -1)
            {
                mxi = x1[i];
                myi = y1[i];
                Mxi = x2[i];
                Myi = y2[i];
            }
            else
            {
                mxi = x2[i];
                myi = y2[i];
                Mxi = x1[i];
                Myi = y1[i];
            }

            if (comparerDoubles(&x1[j], &x2[j]) == -1)
            {
                mxj = x1[j];
                myj = y1[j];
                Mxj = x2[j];
                Myj = y2[j];
            }
            else
            {
                mxj = x2[j];
                myj = y2[j];
                Mxj = x1[j];
                Myj = y1[j];
            }

            if (comparerDoubles(&Mxi, &mxj) == -1  ||  comparerDoubles(&mxi, &Mxj) == 1)
                continue; /* Les 2 segments ne ce chevauchent pas */


            if (comparerDoubles(&mxi, &mxj) == 1)
            {
                grf[0 + 4 * compte] = mxi;
                grf[1 + 4 * compte] = myi;
            }
            else
            {
                grf[0 + 4 * compte] = mxj;
                grf[1 + 4 * compte] = myj;
            }

            if (comparerDoubles(&Mxi, &Mxj) == -1)
            {
                grf[2 + 4 * compte] = Mxi;
                grf[3 + 4 * compte] = Myi;
            }
            else
            {
                grf[2 + 4 * compte] = Mxj;
                grf[3 + 4 * compte] = Myj;
            }
            compte++;
        }

    return;
}

