% notepad
% without sincos
 //<loop> Loop body line 60, nesting depth: 1, estimated iterations: unknown
	.loc	16	98	0
	ld.const.f32 	%f13, [%rd1+20];
	ld.const.f32 	%f14, [%rd1+4];
	ld.const.f32 	%f15, [%rd1+0];
	ld.const.f32 	%f16, [%rd1+8];
	ld.const.f32 	%f17, [%rd1+16];
	sub.ftz.f32 	%f18, %f5, %f14;
	sub.ftz.f32 	%f19, %f8, %f15;
	sub.ftz.f32 	%f20, %f10, %f16;
	mul.ftz.f32 	%f21, %f18, %f18;
	fma.rn.ftz.f32 	%f22, %f19, %f19, %f21;
	fma.rn.ftz.f32 	%f23, %f20, %f20, %f22;
	sqrt.approx.ftz.f32 	%f24, %f23;
	mul.ftz.f32 	%f25, %f9, %f24;
	sub.ftz.f32 	%f26, %f13, %f25;
	mov.f32 	%f27, 0f40c90fdb;    	// 6.28319
	mul.ftz.f32 	%f28, %f24, %f27;
	div.approx.ftz.f32 	%f29, %f17, %f28;
	cos.approx.ftz.f32 	%f30, %f26;
	fma.rn.ftz.f32 	%f12, %f30, %f29, %f12;
	.loc	16	100	0
	sin.approx.ftz.f32 	%f31, %f26;
	fma.rn.ftz.f32 	%f11, %f31, %f29, %f11;
	add.u32 	%r23, %r23, 6;
	add.u64 	%rd1, %rd1, 24;
	setp.gt.u32 	%p5, %r14, %r23;
	@%p5 bra 	$Lt_1_4354;
	bra.uni 	$Lt_1_3842;
    
    
    