.global _start

.hello.str:
    .ascii "12345678\n"

.text

_start:
    movq %rsp, %rbp

    movq $1, %rax
    movq $1, %rdi
    leaq .hello.str, %rsi
    movq $9, %rdx
    syscall

    # exit
    movq $60, %rax
    movq $0, %rdi
    syscall


    pop %rbp