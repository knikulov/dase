# primitive power function
function pow(x, p) {
	res = 1
	for (i = p; i > 0; i--)
		res *= x;

	return res
}

# generates len samples fft code
function gen_fft(n) {
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

	print "	/* conversion to q15 format*/"
	for (i = 0; i < n; i++) {
		print "	X->b[DASE_RE] = (*fp++) * 32768.0;"
		print "	X->b[DASE_IM] = (*fp++) * 32768.0;"
		print "	X++;"
	}

	printf ("\n")
	print "	/* bit reversal sorting */"
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
	print "	/* transform itself */"

	for (frame = 1; frame < n; frame *= 2) {
		arg = 3.14159265358979 / frame;

		print "	/* frame size = " frame " */"
		print "	w.b[DASE_RE] = " 32768.0 * cos(arg) ";"
		print "	w.b[DASE_IM] = inv * " 32768.0 * sin(arg) ";"
		printf ("\n")
		print "	p = X;"
		print "	p2 = X + " frame ";"
		printf ("\n")

		for (i = 0; i < n; i += (2*frame)) {
			fn = i / (2*frame)
			print "	/* frame number = " fn " */"
			print "	s.b[DASE_RE] = 32767.0;"
			print "	s.b[DASE_IM] = 0.0;"
			printf ("\n")

			for (j = 0; j < frame; j++) {
				print "	/* Fourier butterfly */"
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

BEGIN {
	print "#include <stdio.h>"
	print "#include <math.h>"
	printf ("\n")
	print "#include \"dase.h\""
	printf ("\n")
}

NF != 1 {}

(NF == 1) && ($1 > 1) {
	n = pow(2, $1)
	gen_fft(n)
}

