	.file	"stencil.c_compilable.c"
	.section	.text.startup,"ax",@progbits
	.p2align 4,,15
	.globl	main
	.type	main, @function
main:
.LFB5:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movl	$10, %edx
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%r10
	pushq	%rbx
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 10, -56
	.cfi_offset 3, -64
	movq	%rsi, %rbx
	subq	$176, %rsp
	movq	24(%rsi), %rdi
	xorl	%esi, %esi
	call	strtol
	movq	16(%rbx), %rdi
	xorl	%esi, %esi
	movl	$10, %edx
	movq	%rax, %r13
	call	strtol
	movl	%r13d, %r15d
	movq	8(%rbx), %rdi
	xorl	%esi, %esi
	movq	%rax, %r12
	movl	$10, %edx
	call	strtol
	imull	%r12d, %r15d
	movl	$32, %esi
	leaq	-56(%rbp), %rdi
	movq	%rax, %rbx
	imull	%eax, %r15d
	movslq	%r15d, %r14
	salq	$3, %r14
	movq	%r14, %rdx
	call	posix_memalign
	testl	%eax, %eax
	je	.L2
	movq	$0, -56(%rbp)
.L2:
	movq	-56(%rbp), %rax
	movq	%rax, -184(%rbp)
	testl	%r15d, %r15d
	jle	.L15
	shrq	$3, %rax
	negq	%rax
	andl	$3, %eax
	cmpl	%r15d, %eax
	cmova	%r15d, %eax
	cmpl	$6, %r15d
	jg	.L108
	movl	%r15d, %eax
.L6:
	vmovsd	.LC0(%rip), %xmm0
	movq	-184(%rbp), %rsi
	vmovsd	%xmm0, (%rsi)
	cmpl	$1, %eax
	je	.L54
	vmovsd	%xmm0, 8(%rsi)
	cmpl	$2, %eax
	je	.L55
	vmovsd	%xmm0, 16(%rsi)
	cmpl	$3, %eax
	je	.L56
	vmovsd	%xmm0, 24(%rsi)
	cmpl	$4, %eax
	je	.L57
	vmovsd	%xmm0, 32(%rsi)
	cmpl	$6, %eax
	jne	.L58
	vmovsd	%xmm0, 40(%rsi)
	movl	$6, %edx
.L8:
	cmpl	%r15d, %eax
	je	.L15
.L7:
	movl	%r15d, %edi
	leal	-1(%r15), %esi
	movl	%eax, %r9d
	subl	%eax, %edi
	subl	%eax, %esi
	leal	-4(%rdi), %ecx
	shrl	$2, %ecx
	addl	$1, %ecx
	leal	0(,%rcx,4), %r8d
	cmpl	$2, %esi
	jbe	.L10
	vmovapd	.LC1(%rip), %ymm0
	movq	-184(%rbp), %rax
	leaq	(%rax,%r9,8), %rsi
	xorl	%eax, %eax
.L11:
	addl	$1, %eax
	vmovapd	%ymm0, (%rsi)
	addq	$32, %rsi
	cmpl	%eax, %ecx
	ja	.L11
	addl	%r8d, %edx
	cmpl	%edi, %r8d
	je	.L96
	vzeroupper
.L10:
	vmovsd	.LC0(%rip), %xmm0
	movq	-184(%rbp), %rsi
	movslq	%edx, %rax
	vmovsd	%xmm0, (%rsi,%rax,8)
	leal	1(%rdx), %eax
	cmpl	%eax, %r15d
	jle	.L15
	cltq
	vmovsd	%xmm0, (%rsi,%rax,8)
	leal	2(%rdx), %eax
	cmpl	%eax, %r15d
	jle	.L15
	movq	-184(%rbp), %rsi
	cltq
	vmovsd	%xmm0, (%rsi,%rax,8)
.L15:
	movl	var_false(%rip), %ecx
	testl	%ecx, %ecx
	je	.L5
	movq	-184(%rbp), %rdi
	call	dummy
.L5:
	movq	%r14, %rdx
	movl	$32, %esi
	leaq	-56(%rbp), %rdi
	call	posix_memalign
	testl	%eax, %eax
	je	.L16
	movq	$0, -56(%rbp)
.L16:
	movq	-56(%rbp), %rax
	movq	%rax, -192(%rbp)
	testl	%r15d, %r15d
	jle	.L29
	shrq	$3, %rax
	negq	%rax
	andl	$3, %eax
	cmpl	%r15d, %eax
	cmova	%r15d, %eax
	cmpl	$6, %r15d
	jg	.L109
	movl	%r15d, %eax
.L20:
	vmovsd	.LC2(%rip), %xmm0
	movq	-192(%rbp), %rsi
	vmovsd	%xmm0, (%rsi)
	cmpl	$1, %eax
	je	.L61
	vmovsd	%xmm0, 8(%rsi)
	cmpl	$2, %eax
	je	.L62
	vmovsd	%xmm0, 16(%rsi)
	cmpl	$3, %eax
	je	.L63
	vmovsd	%xmm0, 24(%rsi)
	cmpl	$4, %eax
	je	.L64
	vmovsd	%xmm0, 32(%rsi)
	cmpl	$6, %eax
	jne	.L65
	vmovsd	%xmm0, 40(%rsi)
	movl	$6, %edx
.L22:
	cmpl	%r15d, %eax
	je	.L29
.L21:
	movl	%r15d, %edi
	leal	-1(%r15), %esi
	movl	%eax, %r9d
	subl	%eax, %edi
	subl	%eax, %esi
	leal	-4(%rdi), %ecx
	shrl	$2, %ecx
	addl	$1, %ecx
	leal	0(,%rcx,4), %r8d
	cmpl	$2, %esi
	jbe	.L24
	vmovapd	.LC3(%rip), %ymm0
	movq	-192(%rbp), %rax
	leaq	(%rax,%r9,8), %rsi
	xorl	%eax, %eax
.L25:
	addl	$1, %eax
	vmovapd	%ymm0, (%rsi)
	addq	$32, %rsi
	cmpl	%eax, %ecx
	ja	.L25
	addl	%r8d, %edx
	cmpl	%edi, %r8d
	je	.L97
	vzeroupper
.L24:
	vmovsd	.LC2(%rip), %xmm0
	movq	-192(%rbp), %rsi
	movslq	%edx, %rax
	vmovsd	%xmm0, (%rsi,%rax,8)
	leal	1(%rdx), %eax
	cmpl	%eax, %r15d
	jle	.L29
	cltq
	vmovsd	%xmm0, (%rsi,%rax,8)
	leal	2(%rdx), %eax
	cmpl	%eax, %r15d
	jle	.L29
	movq	-192(%rbp), %rsi
	cltq
	vmovsd	%xmm0, (%rsi,%rax,8)
.L29:
	movl	var_false(%rip), %edx
	testl	%edx, %edx
	je	.L19
	movq	-192(%rbp), %rdi
	call	dummy
.L19:
	leal	(%r12,%r12,2), %r14d
	movl	$32, %esi
	leaq	-56(%rbp), %rdi
	imull	%r13d, %r14d
	imull	%ebx, %r14d
	movslq	%r14d, %rdx
	salq	$3, %rdx
	call	posix_memalign
	testl	%eax, %eax
	je	.L30
	movq	$0, -56(%rbp)
.L30:
	movq	-56(%rbp), %rax
	movq	%rax, -200(%rbp)
	testl	%r14d, %r14d
	jle	.L43
	shrq	$3, %rax
	negq	%rax
	andl	$3, %eax
	cmpl	%r14d, %eax
	cmova	%r14d, %eax
	cmpl	$6, %r14d
	jg	.L110
	movl	%r14d, %eax
.L34:
	vmovsd	.LC4(%rip), %xmm0
	movq	-200(%rbp), %rsi
	vmovsd	%xmm0, (%rsi)
	cmpl	$1, %eax
	je	.L68
	vmovsd	%xmm0, 8(%rsi)
	cmpl	$2, %eax
	je	.L69
	vmovsd	%xmm0, 16(%rsi)
	cmpl	$3, %eax
	je	.L70
	vmovsd	%xmm0, 24(%rsi)
	cmpl	$4, %eax
	je	.L71
	vmovsd	%xmm0, 32(%rsi)
	cmpl	$6, %eax
	jne	.L72
	vmovsd	%xmm0, 40(%rsi)
	movl	$6, %edx
.L36:
	cmpl	%r14d, %eax
	je	.L43
.L35:
	movl	%r14d, %edi
	leal	-1(%r14), %r9d
	movl	%eax, %esi
	subl	%eax, %edi
	subl	%eax, %r9d
	leal	-4(%rdi), %ecx
	shrl	$2, %ecx
	addl	$1, %ecx
	leal	0(,%rcx,4), %r8d
	cmpl	$2, %r9d
	jbe	.L38
	vmovapd	.LC5(%rip), %ymm0
	movq	-200(%rbp), %rax
	leaq	(%rax,%rsi,8), %rsi
	xorl	%eax, %eax
.L39:
	addl	$1, %eax
	vmovapd	%ymm0, (%rsi)
	addq	$32, %rsi
	cmpl	%eax, %ecx
	ja	.L39
	addl	%r8d, %edx
	cmpl	%edi, %r8d
	je	.L98
	vzeroupper
.L38:
	vmovsd	.LC4(%rip), %xmm0
	movq	-200(%rbp), %rsi
	movslq	%edx, %rax
	vmovsd	%xmm0, (%rsi,%rax,8)
	leal	1(%rdx), %eax
	cmpl	%eax, %r14d
	jle	.L43
	cltq
	vmovsd	%xmm0, (%rsi,%rax,8)
	leal	2(%rdx), %eax
	cmpl	%eax, %r14d
	jle	.L43
	movq	-200(%rbp), %rsi
	cltq
	vmovsd	%xmm0, (%rsi,%rax,8)
.L43:
	movl	var_false(%rip), %eax
	movl	%eax, -224(%rbp)
	testl	%eax, %eax
	je	.L33
	movq	-200(%rbp), %rdi
	call	dummy
	movl	var_false(%rip), %eax
	movl	%eax, -224(%rbp)
.L33:
	leal	-2(%r12), %eax
	movl	%eax, -176(%rbp)
	cmpl	$2, %eax
	jle	.L44
	leal	-2(%r13), %eax
	movl	%r12d, %esi
	movl	%eax, %edx
	movl	%eax, -96(%rbp)
	leal	-2(%rbx), %eax
	movl	%eax, -212(%rbp)
	movl	%r13d, %eax
	imull	%ebx, %eax
	imull	%eax, %esi
	movl	%eax, -160(%rbp)
	movl	%esi, %edi
	movl	%esi, -220(%rbp)
	movl	%eax, %esi
	addl	%eax, %eax
	movl	%eax, %ecx
	movl	%eax, -216(%rbp)
	cmpl	$2, %edx
	jle	.L44
	movl	%edi, %eax
	movslq	%esi, %r15
	movl	%ecx, -168(%rbp)
	movslq	%ecx, %r14
	addl	%eax, %eax
	movl	%r15d, -164(%rbp)
	movl	%eax, -172(%rbp)
	movslq	%ebx, %rax
	movq	%rax, -112(%rbp)
	salq	$3, %rax
	movq	%rax, -136(%rbp)
	leal	(%rbx,%rbx), %eax
	movl	$2, -156(%rbp)
	movslq	%eax, %rdi
	addl	%ebx, %eax
	cltq
	movq	%rdi, -208(%rbp)
	movq	%rax, -144(%rbp)
	leal	0(,%rbx,4), %eax
	cltq
	movq	%rax, -152(%rbp)
	leal	-5(%rbx), %eax
	addq	$1, %rax
	movq	%rax, -72(%rbp)
.L45:
	addl	$1, -156(%rbp)
	cmpl	$2, -212(%rbp)
	jle	.L47
	movq	-208(%rbp), %rcx
	movslq	-172(%rbp), %rax
	movl	$2, -92(%rbp)
	movq	-200(%rbp), %rsi
	movl	-168(%rbp), %edx
	addq	%rcx, %rax
	movq	-184(%rbp), %rbx
	movslq	-164(%rbp), %r13
	movq	%rcx, -80(%rbp)
	leaq	(%rsi,%rax,8), %rax
	xorl	%esi, %esi
	movq	%rax, -104(%rbp)
	movslq	%edx, %rax
	movq	%r13, %rdi
	leaq	(%rcx,%rax), %r9
	movl	%edi, %r11d
	subl	-160(%rbp), %r11d
	salq	$3, %r9
	movslq	%r11d, %r11
	addq	%r9, %rbx
	subq	%rax, %r11
	addq	-192(%rbp), %r9
	movq	%rbx, -88(%rbp)
	movl	-216(%rbp), %ebx
	movl	%ebx, %r12d
	addl	%edx, %ebx
	addl	%r13d, %r12d
	movslq	%ebx, %rbx
	subq	%rax, %r13
	movslq	%r12d, %r12
	subq	%rax, %rbx
	subq	%rax, %r12
	.p2align 4,,10
	.p2align 3
.L48:
	movq	-144(%rbp), %rax
	movq	-80(%rbp), %r10
	xorl	%ecx, %ecx
	movq	-152(%rbp), %rdi
	addl	$1, -92(%rbp)
	leaq	(%rsi,%rax), %r8
	movq	-88(%rbp), %rax
	movq	-104(%rbp), %rdx
	addq	%rsi, %rdi
	subq	%r10, %r8
	movq	%rax, -120(%rbp)
	movq	-112(%rbp), %rax
	subq	%r10, %rdi
	addq	%rsi, %rax
	subq	-80(%rbp), %rsi
	movq	%rax, -128(%rbp)
	subq	%r10, %rax
	movq	%rax, %r10
	movq	-120(%rbp), %rax
	.p2align 4,,10
	.p2align 3
        movl      $111, %ebx # INSERTED BY KERNCRAFT IACA MARKER UTILITY
        .byte     100        # INSERTED BY KERNCRAFT IACA MARKER UTILITY
        .byte     103        # INSERTED BY KERNCRAFT IACA MARKER UTILITY
        .byte     144        # INSERTED BY KERNCRAFT IACA MARKER UTILITY
.L46:
	vmovsd	16(%rax,%r13,8), %xmm1
	vaddsd	16(%rax,%r12,8), %xmm1, %xmm2
	vmovsd	8(%rax), %xmm1
	vaddsd	24(%rax), %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm1
	vmovsd	16(%rax,%r10,8), %xmm0
	vaddsd	16(%rax,%r8,8), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	16(%rdx), %xmm1
	vmulsd	16(%rdx,%r15,8), %xmm0, %xmm0
	vmulsd	16(%rax), %xmm1, %xmm2
	vaddsd	%xmm2, %xmm0, %xmm1
	vmovsd	16(%rax,%r11,8), %xmm2
	vaddsd	16(%rax,%rbx,8), %xmm2, %xmm3
	vmovsd	(%rax), %xmm2
	vaddsd	32(%rax), %xmm2, %xmm0
	vaddsd	%xmm0, %xmm3, %xmm2
	vmovsd	16(%rax,%rsi,8), %xmm0
	vaddsd	16(%rax,%rdi,8), %xmm0, %xmm0
	addq	$8, %rax
	vaddsd	%xmm0, %xmm2, %xmm0
	vmulsd	16(%rdx,%r14,8), %xmm0, %xmm0
	addq	$8, %rdx
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, 16(%r9,%rcx,8)
	addq	$1, %rcx
	cmpq	-72(%rbp), %rcx
	jne	.L46
        movl      $222, %ebx # INSERTED BY KERNCRAFT IACA MARKER UTILITY
        .byte     100        # INSERTED BY KERNCRAFT IACA MARKER UTILITY
        .byte     103        # INSERTED BY KERNCRAFT IACA MARKER UTILITY
        .byte     144        # INSERTED BY KERNCRAFT IACA MARKER UTILITY
	movq	-136(%rbp), %rax
	movq	-112(%rbp), %rdi
	addq	%rax, -104(%rbp)
	movq	-128(%rbp), %rsi
	addq	%rdi, -80(%rbp)
	addq	%rax, %r9
	movl	-92(%rbp), %edi
	addq	%rax, -88(%rbp)
	cmpl	%edi, -96(%rbp)
	jg	.L48
.L47:
	movl	-160(%rbp), %ebx
	addl	%ebx, -164(%rbp)
	addl	%ebx, -168(%rbp)
	movl	-220(%rbp), %ebx
	addl	%ebx, -172(%rbp)
	movl	-156(%rbp), %ebx
	cmpl	%ebx, -176(%rbp)
	jg	.L45
.L44:
	movl	-224(%rbp), %eax
	testl	%eax, %eax
	jne	.L111
.L99:
	addq	$176, %rsp
	xorl	%eax, %eax
	popq	%rbx
	popq	%r10
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
.L97:
	.cfi_restore_state
	vzeroupper
	jmp	.L29
.L96:
	vzeroupper
	jmp	.L15
.L98:
	vzeroupper
	jmp	.L43
.L111:
	movq	-184(%rbp), %rdi
	call	dummy
	cmpl	$0, var_false(%rip)
	je	.L99
	movq	-192(%rbp), %rdi
	call	dummy
	cmpl	$0, var_false(%rip)
	je	.L99
	movq	-200(%rbp), %rdi
	call	dummy
	jmp	.L99
.L109:
	xorl	%edx, %edx
	testl	%eax, %eax
	je	.L21
	jmp	.L20
.L110:
	xorl	%edx, %edx
	testl	%eax, %eax
	je	.L35
	jmp	.L34
.L108:
	xorl	%edx, %edx
	testl	%eax, %eax
	je	.L7
	jmp	.L6
.L54:
	movl	$1, %edx
	jmp	.L8
.L56:
	movl	$3, %edx
	jmp	.L8
.L55:
	movl	$2, %edx
	jmp	.L8
.L65:
	movl	$5, %edx
	jmp	.L22
.L64:
	movl	$4, %edx
	jmp	.L22
.L58:
	movl	$5, %edx
	jmp	.L8
.L57:
	movl	$4, %edx
	jmp	.L8
.L61:
	movl	$1, %edx
	jmp	.L22
.L68:
	movl	$1, %edx
	jmp	.L36
.L70:
	movl	$3, %edx
	jmp	.L36
.L69:
	movl	$2, %edx
	jmp	.L36
.L72:
	movl	$5, %edx
	jmp	.L36
.L71:
	movl	$4, %edx
	jmp	.L36
.L63:
	movl	$3, %edx
	jmp	.L22
.L62:
	movl	$2, %edx
	jmp	.L22
	.cfi_endproc
.LFE5:
	.size	main, .-main
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC0:
	.long	3961705502
	.long	1071636094
	.section	.rodata.cst32,"aM",@progbits,32
	.align 32
.LC1:
	.long	3961705502
	.long	1071636094
	.long	3961705502
	.long	1071636094
	.long	3961705502
	.long	1071636094
	.long	3961705502
	.long	1071636094
	.section	.rodata.cst8
	.align 8
.LC2:
	.long	424680910
	.long	1071288493
	.section	.rodata.cst32
	.align 32
.LC3:
	.long	424680910
	.long	1071288493
	.long	424680910
	.long	1071288493
	.long	424680910
	.long	1071288493
	.long	424680910
	.long	1071288493
	.section	.rodata.cst8
	.align 8
.LC4:
	.long	3440069995
	.long	1072191488
	.section	.rodata.cst32
	.align 32
.LC5:
	.long	3440069995
	.long	1072191488
	.long	3440069995
	.long	1072191488
	.long	3440069995
	.long	1072191488
	.long	3440069995
	.long	1072191488
	.ident	"GCC: (GNU) 6.3.0"
	.section	.note.GNU-stack,"",@progbits
