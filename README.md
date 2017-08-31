# KASPar_Primer_Deginer

The Perl script is written to automatically design the KASPar Primer for SNP. Firstly, a flanking sequence fasta file must be created using the **Exteact\_Flanking\_Seq.pl** script. The flanking length is suggested to be 300 bp. And then, we can obtain the unique KASPar primer using the script. 

Software requirements:

>**Perl;blastn;primer3_core;**

Usage:

>$ perl KASP.PrimerDesign.Improved.pl -h
>
>  Usage:
>  
>  Run by typing: perl KASP.PrimerDesign.Improved.pl -cfg [Parameter config file] -input [Flanking Seq fasta file] -out [Primer output]
>  
    Required params:
        -c|cfg                                                  [s]     Parameter config file
        -i|input                                                [s]     Flanking Seq fasta file
        -o|out                                                  [s]     Primer output
    Example: perl KASP.PrimerDesign.Improved.pl -cfg [Parameter config file] -input [Flanking Seq fasta file] -out [Primer output]


Output:

![Result](Output.png)
