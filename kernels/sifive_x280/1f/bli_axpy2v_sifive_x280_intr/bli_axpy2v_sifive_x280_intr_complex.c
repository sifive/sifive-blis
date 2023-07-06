/*

   BLIS
   An object-based framework for developing high-performance BLAS-like
   libraries.

   Copyright (C) 2023, SiFive, Inc.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    - Neither the name(s) of the copyright holder(s) nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

// clang-format off
#ifdef AXPY2V

AXPY2V(PRECISION_CHAR, void)
{
    // Computes z := z + alphax * conjx(x) + alphay * conjy(y)
    const DATATYPE* restrict alphax = alphax_;
    const DATATYPE* restrict alphay = alphay_;
    const DATATYPE* restrict x = x_;
    const DATATYPE* restrict y = y_;
    DATATYPE* restrict z = z_;

    if (n <= 0)
        return;

    size_t avl = n;

    while (avl) {
        size_t vl = VSETVL(PREC, LMUL)(avl);
        RVV_TYPE_F_X2(PREC, LMUL) xvec, yvec, zvec;

        if (incx == 1)
            xvec = VLSEG2_V_F(PREC, LMUL)( (BASE_DT*) x, vl);
        else
            xvec = VLSSEG2_V_F(PREC, LMUL)( (BASE_DT*) x, 2*FLT_SIZE*incx, vl);

        if (incy == 1)
            yvec = VLSEG2_V_F(PREC, LMUL)( (BASE_DT*) y, vl);
        else
            yvec = VLSSEG2_V_F(PREC, LMUL)( (BASE_DT*) y, 2*FLT_SIZE*incy, vl);

        if (incz == 1)
            zvec = VLSEG2_V_F(PREC, LMUL)( (BASE_DT*) z, vl);
        else
            zvec = VLSSEG2_V_F(PREC, LMUL)( (BASE_DT*) z, 2*FLT_SIZE*incz, vl);

        RVV_TYPE_F(PREC, LMUL) xvec_real = RVV_GET_REAL(PREC, LMUL, xvec);
        RVV_TYPE_F(PREC, LMUL) xvec_imag = RVV_GET_IMAG(PREC, LMUL, xvec);
        RVV_TYPE_F(PREC, LMUL) yvec_real = RVV_GET_REAL(PREC, LMUL, yvec);
        RVV_TYPE_F(PREC, LMUL) yvec_imag = RVV_GET_IMAG(PREC, LMUL, yvec);
        RVV_TYPE_F(PREC, LMUL) zvec_real = RVV_GET_REAL(PREC, LMUL, zvec);
        RVV_TYPE_F(PREC, LMUL) zvec_imag = RVV_GET_IMAG(PREC, LMUL, zvec);

        //  + alphax * conjx(x)
        zvec_real = VFMACC_VF(PREC, LMUL)( zvec_real, alphax->real, xvec_real, vl);
        zvec_imag = VFMACC_VF(PREC, LMUL)( zvec_imag, alphax->imag, xvec_real, vl);
        if (conjx == BLIS_NO_CONJUGATE){
            zvec_real = VFNMSAC_VF(PREC, LMUL)(zvec_real, alphax->imag, xvec_imag, vl);
            zvec_imag = VFMACC_VF(PREC, LMUL)( zvec_imag, alphax->real, xvec_imag, vl);
        } else {
            zvec_real = VFMACC_VF(PREC, LMUL)( zvec_real, alphax->imag, xvec_imag, vl);
            zvec_imag = VFNMSAC_VF(PREC, LMUL)(zvec_imag, alphax->real, xvec_imag, vl);
        }

        //  + alphay * conjy(y)
        zvec_real = VFMACC_VF(PREC, LMUL)( zvec_real, alphay->real, yvec_real, vl);
        zvec_imag = VFMACC_VF(PREC, LMUL)( zvec_imag, alphay->imag, yvec_real, vl);
        if (conjy == BLIS_NO_CONJUGATE){
            zvec_real = VFNMSAC_VF(PREC, LMUL)(zvec_real, alphay->imag, yvec_imag, vl);
            zvec_imag = VFMACC_VF(PREC, LMUL)( zvec_imag, alphay->real, yvec_imag, vl);
        } else {
            zvec_real = VFMACC_VF(PREC, LMUL)( zvec_real, alphay->imag, yvec_imag, vl);
            zvec_imag = VFNMSAC_VF(PREC, LMUL)(zvec_imag, alphay->real, yvec_imag, vl);
        }

        RVV_SET_REAL(PREC, LMUL, zvec, zvec_real);
        RVV_SET_IMAG(PREC, LMUL, zvec, zvec_imag);

        if (incz == 1)
            VSSEG2_V_F(PREC, LMUL)( (BASE_DT*) z, zvec, vl);
        else
            VSSSEG2_V_F(PREC, LMUL)((BASE_DT*) z, 2*FLT_SIZE*incz, zvec, vl);

        x += vl*incx;
        y += vl*incy;
        z += vl*incz;
        avl -= vl;
    }

}

#endif // AXPY2V
