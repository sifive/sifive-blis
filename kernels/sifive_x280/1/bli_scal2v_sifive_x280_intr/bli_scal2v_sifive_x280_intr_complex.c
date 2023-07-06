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
#ifdef SCAL2V

SCAL2V(PRECISION_CHAR, void)
{
    // Computes y = alpha * conjx(x)
    const DATATYPE* restrict alpha = alpha_;
    const DATATYPE* restrict x = x_;
    DATATYPE* restrict y = y_;

    if (n <= 0) return;
    if (alpha->real == 0 && alpha->imag == 0) {
        SETV(PRECISION_CHAR)(BLIS_NO_CONJUGATE, n, alpha, y, incy, cntx);
        return;
    }

    if (alpha->real == 1 && alpha->imag == 0) {
        COPYV(PRECISION_CHAR)(conjx, n, x, incx, y, incy, cntx);
        return;
    }

    size_t avl = n;
    while (avl) {
        size_t vl = VSETVL(PREC, LMUL)(avl);
        RVV_TYPE_F_X2(PREC, LMUL) xvec, yvec;

        if (incx == 1)
            xvec = VLSEG2_V_F(PREC, LMUL)( (BASE_DT*) x, vl);
        else
            xvec = VLSSEG2_V_F(PREC, LMUL)( (BASE_DT*) x, 2*FLT_SIZE*incx, vl);

        RVV_TYPE_F(PREC, LMUL) xvec_real = RVV_GET_REAL(PREC, LMUL, xvec);
        RVV_TYPE_F(PREC, LMUL) xvec_imag = RVV_GET_IMAG(PREC, LMUL, xvec);
        RVV_TYPE_F(PREC, LMUL) yvec_real = RVV_GET_REAL(PREC, LMUL, yvec);
        RVV_TYPE_F(PREC, LMUL) yvec_imag = RVV_GET_IMAG(PREC, LMUL, yvec);

        yvec_real = VFMUL_VF(PREC, LMUL)(xvec_real, alpha->real, vl);
        yvec_imag = VFMUL_VF(PREC, LMUL)(xvec_real, alpha->imag, vl);
        if (conjx == BLIS_NO_CONJUGATE) {
            yvec_real = VFNMSAC_VF(PREC, LMUL)(yvec_real, alpha->imag, xvec_imag, vl);
            yvec_imag = VFMACC_VF( PREC, LMUL)(yvec_imag, alpha->real, xvec_imag, vl);
        } else {
            yvec_real = VFMACC_VF( PREC, LMUL)(yvec_real, alpha->imag, xvec_imag, vl);
            yvec_imag = VFNMSAC_VF(PREC, LMUL)(yvec_imag, alpha->real, xvec_imag, vl);
        }

        RVV_SET_REAL(PREC, LMUL, yvec, yvec_real);
        RVV_SET_IMAG(PREC, LMUL, yvec, yvec_imag);
#pragma GCC diagnostic ignored "-Wuninitialized"

        if (incy == 1)
            VSSEG2_V_F(PREC, LMUL)( (BASE_DT*) y, yvec, vl);
        else
            VSSSEG2_V_F(PREC, LMUL)((BASE_DT*) y, 2*FLT_SIZE*incy, yvec, vl);

        x += vl*incx;
        y += vl*incy;
        avl -= vl;
    }

}

#endif // SCAL2V
