	.file	"dummy.c"
	.text
	.p2align 4,,15
	.globl	dummy
	.type	dummy, @function
dummy:
.LFB0:
	.cfi_startproc
	ret
	.cfi_endproc
.LFE0:
	.size	dummy, .-dummy
	.globl	var_false
	.bss
	.align 4
	.type	var_false, @object
	.size	var_false, 4
var_false:
	.zero	4
	.ident	"GCC: (GNU) 6.3.0"
	.section	.note.GNU-stack,"",@progbits
