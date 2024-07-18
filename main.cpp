#include <stdlib.h>

#include "automorphisms.h"

int main() {
    MODULE* module = new_module_info(N, FFT64);

    Int64VecN skey_raw;
    SVP_PPOL* skey = new_svp_ppol(module);

    Int64VecN mu[message_limbs];

    int64_t p = 3;

    VMP_PMAT* ks_a = new_vmp_pmat(module, autom_nrows, autom_ncols);
    VMP_PMAT* ks_b = new_vmp_pmat(module, autom_nrows, autom_ncols);

    Int64VecN a[ell];
    Int64VecN b[ell];

    Int64VecN autom_a[ell];
    Int64VecN autom_b[ell];

    // generate a secret key
    random_binary(N, skey_raw);
    svp_prepare(module, skey, skey_raw);

    // generate the autom ks
    create_keyswitch(module, p, K,  //
                     ks_a, ks_b,    //
                     skey_raw, skey);

    // generate a random message
    for (uint64_t i = 0; i < message_limbs; ++i) {
        random_centered_reduced(N, K, mu[i]);
    }

    // encrypt it (noise level = K*ell)
    rlwe_encrypt(module, K, *a, *b, ell, *mu, message_limbs, skey);

    auto start = std::chrono::steady_clock::now();

    for (int i = 0; i < 50000; i++) {
        apply_automorphism(module, p, K,        //
                           *autom_a, *autom_b,  //
                           *a, *b,              //
                           ks_a, ks_b);         //
    }

    auto finish = std::chrono::steady_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::duration<double, std::micro>>(finish - start)
            .count();
    std::cout << elapsed / 50000 << "us" << std::endl;

    delete_svp_ppol(skey);
    delete_vmp_pmat(ks_b);
    delete_vmp_pmat(ks_a);
    delete_module_info(module);
}