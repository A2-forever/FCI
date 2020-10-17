#ifndef CI_CPP
#define CI_CPP
#include "FCI.h"

double CI::H(Slater_det &k1, Slater_det &k2)
{
    vector<int> Num(2);//位置不同的电子的数量
    vector<int> index(2 * 2 * nOrb);//2个组态，2种自旋，最多nOrb个轨道

    find(k1, k2, Num, index);
    
    double h1 = 0; // one-electron contribution
    double h2 = 0; // two-e+lectron contribution
    
    int n_excit_a = Num[0];
    int n_excit_b = Num[1];
    if(n_excit_a + n_excit_b <= 2){
        // diagonal element
        if(n_excit_a == 0 && n_excit_b == 0){
            // one-electron contribution
            // h1 = [h1e(p, p) for p in occ_list_a] 
            //    + [h1e(p, p) for p in occ_list_b]
            // p in occ_list_a
            for(int p = 0; p < nOrb; p++){
                h1 += k1.Orb(p, 0) * h[p * dim + p];
            }
            // p in occ_list_b
            for(int p = 0; p < nOrb; p++){
                h1 += k1.Orb(p, 1) * h[p * dim + p];
            }
            // two-electron contribution
            // h2 = 0.5 * {[h2e(p, p, r, r) - h2e(p, r, r, p)] for [(p, r in occ_list_a) + (p, r in occ_list_b)+-
            //           + [h2e(p, p, r, r)                    for [(p in occ_list_a, r in occ_list_b) + 
            //                                                      (p in occ_list_b, r in occ_list_a)]}
            // p in occ_list_a
            for(int p = 0; p < nOrb; p++)
            {
                // r in occ_list_a
                for(int r = 0; r < nOrb; r++)
                {
                    h2 += k1.Orb(p, 0) * k1.Orb(r, 0) * (g[p * dim3 + p * dim2 + r * dim + r] - g[p * dim3 + r * dim2 + r * dim + p]);
                }
                // r in occ_list_b
                for(int r = 0; r < nOrb; r++)
                {
                    h2 += k1.Orb(p, 0) * k1.Orb(r, 1) * g[p * dim3 + p * dim2 + r * dim + r];
                }
            }
            // p in occ_list_b
            for(int p = 0; p < nOrb; p++)
            {
                // r in occ_list_a
                for(int r = 0; r < nOrb; r++)
                {
                    h2 += k1.Orb(p, 1) * k1.Orb(r, 0) * g[p * dim3 + p * dim2 + r * dim + r];
                }
                // r in occ_list_b
                for(int r = 0; r < nOrb; r++)
                {
                    h2 += k1.Orb(p, 1) * k1.Orb(r, 1) * (g[p * dim3 + p * dim2 + r * dim + r] - g[p * dim3 + r * dim2 + r * dim + p]);
                }
            }
            h2 *= 0.5;
        }
        // a -> a
        else if(n_excit_a == 1 && n_excit_b == 0){
            int p, q;
            p = index[1 * 2 * dim + 0 * dim + 0];
            q = index[0 * 2 * dim + 0 * dim + 0];

            double sign = k2.gamma(p, 0) * k1.gamma(q, 0);
            // one-electron contribution
            // h1 = sign(p, q) * h1e(p, q)
            h1 = sign * h[p * dim + q];
            
            // two-electron contribution
            // h2 = sign(p, q) * {[h2e(p, q, r, r) - h2e(p, r, r, q)] for r in occ_list_a
            //                   + h2e(p, q, r, r) for r in occ_list_b}
            // r in occ_list_a
            for(int r = 0; r < nOrb; r++)
            {
                h2 += k1.Orb(r, 0) * (g[p * dim3 + q * dim2 + r * dim + r] - g[p * dim3 + r * dim2 + r * dim + q]);
            }            
            // r in occ_list_b
            for(int r = 0; r < nOrb; r++)
            {
                h2 += k1.Orb(r, 1) * g[p * dim3 + q * dim2 + r * dim + r];
            }
            h2 *= sign;
        }
        // b -> b
        else if(n_excit_a == 0 && n_excit_b == 1)
        {
            int p, q;
            p = index[1 * 2 * dim + 1 * dim + 0];
            q = index[0 * 2 * dim + 1 * dim + 0];

            double sign = k2.gamma(p, 1) * k1.gamma(q, 1);
            // one-electron contribution
            // h1 = sign(p, q) * h1e(p, q)
            h1 = sign * h[p * dim + q];
            
            // two-electron contribution
            // h2 = sign(p, q) * {[h2e(p, q, r, r) - h2e(p, r, r, q)] for r in occ_list_b
            //                   + h2e(p, q, r, r) for r in occ_list_a}
            // r in occ_list_b
            for(int r = 0; r < nOrb; r++)
            {
                h2 += k1.Orb(r, 1) * (g[p * dim3 + q * dim2 + r * dim + r] - g[p * dim3 + r * dim2 + r * dim + q]);
            }
            // r in occ_list_a
            for(int r = 0; r < nOrb; r++)
            {
                h2 += k1.Orb(r, 0) * g[p * dim3 + q * dim2 + r * dim + r];
            }
            h2 *= sign;
        }
        // a, a -> a, a
        else if(n_excit_a == 2){
            int p = index[1 * 2 * dim + 0 * dim + 0];
            int r = index[0 * 2 * dim + 0 * dim + 0];
            int q = index[1 * 2 * dim + 0 * dim + 1];
            int s = index[0 * 2 * dim + 0 * dim + 1];

            double sign = k2.gamma(p, 0) * k1.gamma(r, 0) * k2.gamma(q, 0) * k1.gamma(s, 0);
            // two-electron contribution
            // h2 = sign(p, q) * [h2e(p, q, r, s) - h2e(p, s, r, q)]
            h2 = sign * (g[p * dim3 + q * dim2 + r * dim + r] - g[p * dim3 + s * dim2 + r * dim + q]);

           
        }
        // b, b -> b, b
        else if(n_excit_b == 2){
            int p = index[1 * 2 * dim + 1 * dim + 0];
            int r = index[0 * 2 * dim + 1 * dim + 0];
            int q = index[1 * 2 * dim + 1 * dim + 1];
            int s = index[0 * 2 * dim + 1 * dim + 1];
            // make sure p < r and q < s
            //std::cout << "p = " << p << " r = " << r << " q = " << q << " s = " << s << '\n';
            //std::cout << p << " " << q << " " << r << " " << s << " " << _Integrals.h2e(p, q, r, s) << '\n';
            //std::cout << p << " " << s << " " << r << " " << q << " " << _Integrals.h2e(p, s, r, q) << '\n';
            double sign = k2.gamma(p, 1) * k1.gamma(r, 1) * k2.gamma(q, 1) * k1.gamma(s, 1);
            // two-electron contribution
            // h2 = sign(p, q) * [h2e(p, q, r, s) - h2e(p, s, r, q)]
            h2 = sign * (g[p * dim3 + q * dim2 + r * dim + r] - g[p * dim3 + s * dim2 + r * dim + q]);

        }
        // a, b -> a, b
        else if(n_excit_a == 1 && n_excit_b == 1){
            int p = index[1 * 2 * dim + 0 * dim + 0];
            int q = index[0 * 2 * dim + 0 * dim + 0];
            int r = index[1 * 2 * dim + 1 * dim + 0];
            int s = index[0 * 2 * dim + 1 * dim + 0];

            double sign = k2.Orb(p, 0) * k1.Orb(q, 0);
            sign *= k2.Orb(r, 0) * k1.Orb(s, 0);
            // two-electron contribution
            // h2 = sign(p, q) * h2e(p, q, r, s)
            h2 = sign * g[p * dim3 + q * dim2 + r * dim + s];
        }
    }
    return h1 + h2;
}


#endif