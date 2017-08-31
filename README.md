# KASPar_Primer_Deginer

The Perl script is written to automatically design the KASPar (Kompetitive Allele Specific PCR) Primer for SNP. Firstly, a flanking sequence fasta file must be created using the **Exteact\_Flanking\_Seq.pl** script. The flanking length is suggested to be 300 bp. And then, we can obtain the unique KASPar primer using the script. 

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

**Method:**

The KASPar primers consist of a allele-specific forward primer and a common reverse primer. We assumes that the forward primer is a primer with constant length (the maximum length that allowed in the configure file). The reverse primer is designed using the software Primer3. And the forward and reverse primer were botn aligned to the reference genome by use of blast+, and only the unique pairs (at least one is unique) were picked in the last. The suggestive parameter for Primer3 and blast can be found in the configure file here.