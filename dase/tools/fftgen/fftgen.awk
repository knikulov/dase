# generates len samples radix 2 fft code
function gen_fft_radix2(n) {
	print "void"
	print "dase_fft" n "(float *x, float inv)"
	print "{"
	print "	a64 ac0;"
	print "	v2q15_union w, s, sinv, tmp;"
	print "	v2q15_union *X, *p, *p2;"
	print "	float *fp;"
	printf ("\n")

	print "	fp = x;"
	print "	X = (v2q15_union *)x;"
	printf ("\n")

	for (i = 0; i < n; i++) {
		print "	X->b[DASE_RE] = (*fp++) * 32768.0;"
		print "	X->b[DASE_IM] = (*fp++) * 32768.0;"
		print "	X++;"
	}

	printf ("\n")
	print "	X = (v2q15_union *)x;"
	printf ("\n")

	for (i = 0; i < n; i++) {
		s = 0
		for (c = 1; c < n; c *= 2) {
			s *= 2
			if ((int(i / c) % 2) == 1) {
				s++
			}
		}

		if (s > i) {
			print "	tmp = X[" i "];"
			print "	X[" i "] = X[" s "];"
			print "	X[" s "] = tmp;"
		}
	}

	printf ("\n")

	for (frame = 1; frame < n; frame *= 2) {
		arg = 3.14159265358979 / frame;

		print "	w.b[DASE_RE] = " 32768.0 * cos(arg) ";"
		print "	w.b[DASE_IM] = inv * " 32768.0 * sin(arg) ";"
		printf ("\n")
		print "	p = X;"
		print "	p2 = X + " frame ";"
		printf ("\n")

		for (i = 0; i < n; i += (2*frame)) {
			fn = i / (2*frame)
			print "	s.b[DASE_RE] = 32767.0;"
			print "	s.b[DASE_IM] = 0.0;"
			printf ("\n")

			for (j = 0; j < frame; j++) {
				print "	sinv.b[DASE_RE] = s.b[DASE_IM];"
				print "	sinv.b[DASE_IM] = s.b[DASE_RE];"
				print "	ac0 = 0;"
				print "	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);"
				print "	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);"
				print "	ac0 = 0;"
				print "	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);"
				print "	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);"
				print "	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);"
				print "	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);"
				print "	ac0 = 0;"
				print "	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);"
				print "	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);"
				print "	ac0 = 0;"
				print "	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);"
				print "	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);"
				printf ("\n") 
				print "	p++;"
				print "	p2++;"
				printf ("\n")
			}
			print "	p  += " frame ";"
			print "	p2 += " frame ";"
			printf ("\n")
		}
	}

	print "	fp = x + 2*" n ";"
	print "	X = (v2q15_union *)x + " n ";"
	printf ("\n")

	for (i = 0; i < n; i++) {
		print "	--X;"
		print "	*--fp = X->b[DASE_IM] / 32768.0;"
		print "	*--fp = X->b[DASE_RE] / 32768.0;"
	}

	print "}"
	printf ("\n")
}

# generates len samples radix 4 fft code
function gen_fft_radix4(n) {
	gen_fft_radix4_dir(n)
	gen_fft_radix4_inv(n)

	print "void"
	print "dase_fft" n "(float *x, float inv)"
	print "{"
	print "	if (inv == DASE_FFT_DIR)"
	print "		dase_fft" n "_dir(x);"
	print "	else"
	print "		dase_fft" n "_inv(x);"
	print "}"
	printf ("\n")
}

function gen_fft_radix4_dir(n){
	print "void"
	print "dase_fft" n "_dir(float *x)"
	print "{"
	print "	a64 ac0;"
	print "	v2q15_union w0, w1, w2;"
	print "	v2q15_union w0i, w1i, w2i;"
	print "	v2q15_union t0, t1, t2, t3;"
	print "	v2q15_union b0, b1;"
	print "	v2q15_union tmp;"
	print "	v2q15_union *X0, *X1, *X2, *X3;"
	print "	v2q15_union *X;"
	print "	float *fp;"
	printf ("\n")

	print "	fp = x;"
	print "	X = (v2q15_union *)x;"
	printf ("\n")

	for (i = 0; i < n; i++) {
		print "	X->b[DASE_RE] = (*fp++) * 32768.0;"
		print "	X->b[DASE_IM] = (*fp++) * 32768.0;"
		print "	X++;"
	}

	printf ("\n")
	print "	X = (v2q15_union *)x;"
	printf ("\n")

	n2 = n;
	for (k = 1; k < n; k *= 4) {
		n1 = n2
		n2 = n2 / 4
		e = 2 * 3.14159265358979 / n1

		for (j = 0; j < n2; j++) {
			arg0 = j*e;
			arg1 = arg0 + arg0
			arg2 = arg0 + arg1

			print "	w0.b[DASE_RE] = w0i.b[DASE_IM] = " 32767.0 * cos(arg0) ";"
			print "	w0.b[DASE_IM] = w0i.b[DASE_RE] = " 32767.0 * sin(arg0) ";"
			print "	w1.b[DASE_RE] = w1i.b[DASE_IM] = " 32767.0 * cos(arg0) ";"
			print "	w1.b[DASE_IM] = w1i.b[DASE_RE] = " 32767.0 * sin(arg0) ";"
			print "	w2.b[DASE_RE] = w2i.b[DASE_IM] = " 32767.0 * cos(arg0) ";"
			print "	w2.b[DASE_IM] = w2i.b[DASE_RE] = " 32767.0 * sin(arg0) ";"
			printf ("\n")

			for (i = j; i < n; i += n1) {
				i1 = i + n2
				i2 = i1 + n2
				i3 = i2 + n2

				print "	X0 = X + " i ";"
				print "	X1 = X + " i1 ";"
				print "	X2 = X + " i2 ";"
				print "	X3 = X + " i3 ";"
				printf ("\n")
				print "	t0.a = __builtin_mips_addq_s_ph(X0->a, X2->a);"
				print "	t1.a = __builtin_mips_addq_s_ph(X1->a, X3->a);"
				print "	t2.a = __builtin_mips_subq_s_ph(X0->a, X2->a);"
				print "	t3.a = __builtin_mips_subq_s_ph(X1->a, X3->a);"
				printf ("\n")
				print "	b0.a = t2.a;"
				print "	b1.b[DASE_RE] = -t3.b[DASE_IM];"
				print "	b1.b[DASE_IM] =  t3.b[DASE_RE];"
				printf ("\n")
				print "	X0->a = __builtin_mips_addq_s_ph(t0.a, t1.a);"
				print "	t1.a = __builtin_mips_subq_s_ph(t0.a, t1.a);"
				print "	t0.a = __builtin_mips_addq_s_ph(b0.a, b1.a);"
				print "	t2.a = __builtin_mips_subq_s_ph(b0.a, b1.a);"
				printf ("\n")
				print "	ac0 = 0;"
				print "	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, t2.a, w0.a);"
				print "	X1->b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);"
				printf ("\n")
				print "	ac0 = 0;"
				print "	ac0 -= __builtin_mips_mulsaq_s_w_ph(ac0, t2.a, w0i.a);"
				print "	X1->b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);"
				printf ("\n")
				print "	ac0 = 0;"
				print "	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, t1.a, w1.a);"
				print "	X2->b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);"
				printf ("\n")
				print "	ac0 = 0;"
				print "	ac0 -= __builtin_mips_mulsaq_s_w_ph(ac0, t1.a, w1i.a);"
				print "	X2->b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);"
				printf ("\n")
				print "	ac0 = 0;"
				print "	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, t0.a, w2.a);"
				print "	X3->b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);"
				printf ("\n")
				print "	ac0 = 0;"
				print "	ac0 -= __builtin_mips_mulsaq_s_w_ph(ac0, t0.a, w2i.a);"
				print "	X3->b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);"
				printf ("\n")
			}
		}
	}

	print "	X = (v2q15_union *)x;"
	j = 0
	n2 = n / 4
	for (i = 0; i < n-1; i++) {
		if (i < j) {
			print "	tmp = X[" i "];"
			print "	X[" i "] = X[" j "];"
			print "	X[" j "] = tmp;"
		}
		n1 = n2;
		while (j >= 3*n1) {
			j -= 3*n1
			n1 /= 4
		}
		j += n1
	}

	print "}"
	printf ("\n")
}

function gen_fft_radix4_inv(n){
	print "void"
	print "dase_fft" n "_inv(float *x)"
	print "{"
	print "	a64 ac0;"
	print "	v2q15_union w0, w1, w2;"
	print "	v2q15_union w0i, w1i, w2i;"
	print "	v2q15_union t0, t1, t2, t3;"
	print "	v2q15_union b0, b1;"
	print "	v2q15_union tmp;"
	print "	v2q15_union *X0, *X1, *X2, *X3;"
	print "	v2q15_union *X;"
	print "	float *fp;"
	printf ("\n")

	print "	fp = x;"
	print "	X = (v2q15_union *)x;"
	printf ("\n")

	for (i = 0; i < n; i++) {
		print "	X->b[DASE_RE] = (*fp++) * 32768.0;"
		print "	X->b[DASE_IM] = (*fp++) * 32768.0;"
		print "	X++;"
	}

	printf ("\n")
	print "	X = (v2q15_union *)x;"
	printf ("\n")

	n2 = n;
	for (k = 1; k < n; k *= 4) {
		n1 = n2
		n2 = n2 / 4
		e = 2 * 3.14159265358979 / n1

		for (j = 0; j < n2; j++) {
			arg0 = j*e;
			arg1 = arg0 + arg0
			arg2 = arg0 + arg1

			print "	w0.b[DASE_RE] = w0i.b[DASE_IM] = " 32767.0 * cos(arg0) ";"
			print "	w0.b[DASE_IM] = w0i.b[DASE_RE] = " 32767.0 * sin(arg0) ";"
			print "	w1.b[DASE_RE] = w1i.b[DASE_IM] = " 32767.0 * cos(arg0) ";"
			print "	w1.b[DASE_IM] = w1i.b[DASE_RE] = " 32767.0 * sin(arg0) ";"
			print "	w2.b[DASE_RE] = w2i.b[DASE_IM] = " 32767.0 * cos(arg0) ";"
			print "	w2.b[DASE_IM] = w2i.b[DASE_RE] = " 32767.0 * sin(arg0) ";"
			printf ("\n")

			for (i = j; i < n; i += n1) {
				i1 = i + n2
				i2 = i1 + n2
				i3 = i2 + n2

				print "	X0 = X + " i ";"
				print "	X1 = X + " i1 ";"
				print "	X2 = X + " i2 ";"
				print "	X3 = X + " i3 ";"
				printf ("\n")
				print "	t0.a = __builtin_mips_addq_s_ph(X0->a, X2->a);"
				print "	t1.a = __builtin_mips_addq_s_ph(X1->a, X3->a);"
				print "	t2.a = __builtin_mips_subq_s_ph(X0->a, X2->a);"
				print "	t3.a = __builtin_mips_subq_s_ph(X1->a, X3->a);"
				printf ("\n")
				print "	b0.a = t2.a;"
				print "	b1.b[DASE_RE] = -t3.b[DASE_IM];"
				print "	b1.b[DASE_IM] =  t3.b[DASE_RE];"
				printf ("\n")
				print "	X0->a = __builtin_mips_addq_s_ph(t0.a, t1.a);"
				print "	t1.a = __builtin_mips_subq_s_ph(t0.a, t1.a);"
				print "	t0.a = __builtin_mips_addq_s_ph(b0.a, b1.a);"
				print "	t2.a = __builtin_mips_subq_s_ph(b0.a, b1.a);"
				printf ("\n")
				print "	ac0 = 0;"
				print "	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, t2.a, w0.a);"
				print "	X3->b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);"
				printf ("\n")
				print "	ac0 = 0;"
				print "	ac0 -= __builtin_mips_dpaq_s_w_ph(ac0, t2.a, w0i.a);"
				print "	X3->b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);"
				printf ("\n")
				print "	ac0 = 0;"
				print "	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, t1.a, w1.a);"
				print "	X2->b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);"
				printf ("\n")
				print "	ac0 = 0;"
				print "	ac0 -= __builtin_mips_dpaq_s_w_ph(ac0, t1.a, w1i.a);"
				print "	X2->b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);"
				printf ("\n")
				print "	ac0 = 0;"
				print "	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, t0.a, w2.a);"
				print "	X1->b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);"
				printf ("\n")
				print "	ac0 = 0;"
				print "	ac0 -= __builtin_mips_dpaq_s_w_ph(ac0, t0.a, w2i.a);"
				print "	X1->b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);"
				printf ("\n")
			}
		}
	}

	print "	X = (v2q15_union *)x;"
	j = 0
	n2 = n / 4
	for (i = 0; i < n-1; i++) {
		if (i < j) {
			print "	tmp = X[" i "];"
			print "	X[" i "] = X[" j "];"
			print "	X[" j "] = tmp;"
		}
		n1 = n2;
		while (j >= 3*n1) {
			j -= 3*n1
			n1 /= 4
		}
		j += n1
	}

	print "}"
	printf ("\n")
}

BEGIN {
	print "#include \"dase.h\""
	printf ("\n")
}

NF != 1 {}

(NF == 1) && ($1 > 1) {
	if ($1 % 2 == 1) {
		gen_fft_radix2(2^$1)
	} else {
		gen_fft_radix4(2^$1)
	}
}

