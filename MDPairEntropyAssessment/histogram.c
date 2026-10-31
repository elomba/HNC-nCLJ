/* Unordered minimum-image pair counts; input positions must be wrapped. */
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <omp.h>

int pair_histogram(const double *xyz, const int32_t *types, int n,
                   const double *box, const int32_t *pairmap,
                   int nspecies, int nbins, double dr, int nthreads,
                   uint32_t *out) {
    if (n < 2 || nbins < 1 || dr <= 0 || nthreads < 1) return -1;
    int npairs = nspecies * (nspecies + 1) / 2;
    size_t size = (size_t)npairs * nbins;
    uint32_t *local = calloc((size_t)nthreads * size, sizeof(uint32_t));
    if (!local) return -2;
    double limit = nbins * dr, limit2 = limit * limit;
    for (size_t k=0; k<size; ++k) out[k] = 0;
    #pragma omp parallel num_threads(nthreads)
    {
        uint32_t *h = local + (size_t)omp_get_thread_num() * size;
        #pragma omp for schedule(dynamic,16)
        for (int i=0; i<n-1; ++i) {
            double xi=xyz[3*i], yi=xyz[3*i+1], zi=xyz[3*i+2];
            for (int j=i+1; j<n; ++j) {
                double dx=xi-xyz[3*j], dy=yi-xyz[3*j+1], dz=zi-xyz[3*j+2];
                if(dx>box[0]*.5) dx-=box[0]; else if(dx < -box[0]*.5) dx+=box[0];
                if(dy>box[1]*.5) dy-=box[1]; else if(dy < -box[1]*.5) dy+=box[1];
                if(dz>box[2]*.5) dz-=box[2]; else if(dz < -box[2]*.5) dz+=box[2];
                double r2=dx*dx+dy*dy+dz*dz;
                if(r2 < limit2) {
                    int b=(int)(sqrt(r2)/dr);
                    if(b<nbins) ++h[(size_t)pairmap[types[i]*nspecies+types[j]]*nbins+b];
                }
            }
        }
    }
    for(int t=0;t<nthreads;++t)
        for(size_t k=0;k<size;++k) out[k]+=local[(size_t)t*size+k];
    free(local);
    return 0;
}
