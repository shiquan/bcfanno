Interpret the annotations.
==========================

Interpret human genetic variants
--------------------------------

+ Interpret coding variants

    The American College of Medical Genetics and Genomics (ACMG) supply a decision-tree roadmap and recommend using 28 criteria to help the clinical researcher to interpret genetic variants, however no computation approach could interpret the genetic variants directly, that's because gathering information for all the criteria is quite complicated and no specific algorithms for implementing this guidelines specified. See details, please refer to  **ACMG guideline** (https://www.acmg.net/docs/standards_guidelines_for_the_interpretation_of_sequence_variants.pdf).

+ Interpret noncoding/regulatory variants
      
    

+ VarType and HGVSnom


**Variant types**

For genetic variants located in the transcript regions, bcfanno will interpret the amino acid changes and the variant types (**VarType**). One major challange encountered with interpreting VarType is multi transcripts may be encoded by one gene, several VarTypes will be exported and separated with "|" in such a record. Beside VarType, any other transcription-related tags (records), like **HGVSnom**, will also be generated in such a format, which multi transcription-related information seperated with "|" and exported in order.



Supported VarTypes are summaried below.

* *Intron_* :  genetic variant in intron;

.. _Intron:http://www.sequenceontology.org/miso/current_release/term/SO:0001627

* *Noncoding_* :  genetic variant in the exonic region of a noncoding transcription;

.. _Noncoding:http://www.sequenceontology.org/miso/current_release/term/SO:0001792

* *Utr5_* :  genetic variant in the 5' untranslated region (except splice sites);

.. _Utr5:http://www.sequenceontology.org/miso/current_release/term/SO:0001623

* *Utr3_* :  genetic variant in the 3' untranslated region (except splice sites);

.. _Utr3:http://www.sequenceontology.org/miso/current_release/term/SO:0001624

* *Synonymous_* :  genetic variants in coding sequence, resulting in no change to encoded amino acid;

.. _Synonymous:http://www.sequenceontology.org/miso/current_release/term/SO:0001819

* *Missense_* :  genetic variants in coding sequence, resulting in a different amino acid;

.. _Missense:http://www.sequenceontology.org/miso/current_release/term/SO:0001583

* *Nonsense_* :  genetic variants in coding sequence, resulting in a premature stop codon, leading to a shortened transcript, also known as *stop-gain*;

.. _Nonsense:http://www.sequenceontology.org/miso/current_release/term/SO:0001587

* *InframeInsertion_* :  insert bases in the coding sequence, resulting in no reading frame shift;

.. _InframeInsertion:http://www.sequenceontology.org/miso/current_release/term/SO:0001821

* *InframeDeletion_* :  delete bases in the coding sequence, resulting in no reading frame shift;

.. _InframeDeletion:http://www.sequenceontology.org/miso/current_release/term/SO:0001822

* *InframeDelins_* :  delete and insert bases in the coding sequence, resulting in no reading frame shift;

.. _InframeDelins:http://www.sequenceontology.org/miso/current_release/term/SO:0001820

* *Frameshift_* :  delete or insert bases in the coding sequence, resulting in reading frame shift, leading to a incomplete transcript;

.. _Frameshift:http://www.sequenceontology.org/miso/current_release/term/SO:0001589

* *StopLost_* :  genetic variants in the stop codon, resulting in an elongated transcript;

.. _StopLost:http://www.sequenceontology.org/miso/current_release/term/SO:0001578

* *StopRetained_* :  genetic variants in the stop codon but terminator remains;

.. _StopRetained:http://www.sequenceontology.org/miso/current_release/term/SO:0001567

* *SpliceSite_* :  genetic variants in the region of splice site, either within 1-3 bases of the exon or 3-8 bases of the intron, please notice that 1-2 bases of the intron will be interpret as *SpliceDonor* or *SpliceAcceptor*; 

.. _SpliceSite:http://www.sequenceontology.org/miso/current_release/term/SO:0001630

* *SpliceDonor_* :  genetic variants in the 2 base region at the 5' end of an intron;

.. _SpliceDonor:http://www.sequenceontology.org/miso/current_release/term/SO:0001575

* *SpliceAcceptor_* : genetic variants in the 2 base region at the 3' end of an intron;

.. _SpliceAcceptor:http://www.sequenceontology.org/miso/current_release/term/SO:0001574


* *Complex*,  large variants overlapped coding and noncoding regions, or influence more than one gene;

* *NoCall*,  genetic variants in coding region, and the alternative allele is the same with the base in the transcript sequence, the *NoCall* type come from the inconsistance between genome reference and transcription reference, and usually come with high allele frequency in population database;

* *Unknown*, genetic variants in coding region but could not be interpreted by bcfanno, usually account for program bugs.

For any other variants not annotated with bcfanno (empty or no VarType tag) could be interpret as intergenic type.


### Impact order
bcfanno predict variants in order of severity.

* *Complex*

* *Nonsense*, *StopLost*, *StopRetained*, *SpliceDonor*, *SpliceAcceptor*, *SpliceSite*, *Frameshift*

* *InframeDelins*, *InframeDeletion*, *InframeInsertion*, *Missense*

* *Synonymous*

* *Utr5*, *Utr3*, *NoCall*, *Intron*, *Noncoding*

* *Unknown*

## HGVS nomenclature

HGVS is short for Human Genome Variation Society. Nowadays, HGVS nomenclature (**HGVSnom**) is recommended to report and describe sequence variants found in DNA, RNA and protein. bcfanno generate three kinds of HGVS tags to describe the variants in gene regions, **HGVSnom** is the standard HGVS nomenclature, **Oldnom** descibe the variant location in gene without count the UTR regions, **IVSnom** describe the variants in intron. **Oldnom** and **IVSnom** only used to check the record published several years ago. It is recommend to descibe and publish genetic variants in standard format.

Here is a demo:

::
   
   #CHROM	START	END	TGT	HGVSnom	Oldnom	IVSnom
   chr8	37821852	37821853	T/C	NM_000025.2:c.1206-96A>G	NM_000025.2:n.1403-96A>G	NM_000025.2:c.IVS1-96T>C
   chr1	11906067	11906068	A/A	NR_037806.1:n.1479+245A>G|NM_006172.3:c.454T>C(p.Stop152Arg/p.X152R)	NM_006172.3:n.553T>C	NR_037806.1:c.IVS3+245A>G

   
Please notice that one gene may encode more than one transcript, bcfanno will annotate all the transcripts in the refgene databases if no transcript list specified. And the HGVS names of different transcripts will seperated with "|". If two or more transcripts found, the **Gene** tag, **Transcripts** tag and **HGVSnom** tag will generated in same order. However, **Oldnom** will be exported only if coding variants in the transcript and **IVSnom** for intron variants.



Reference
---------

* http://varnomen.hgvs.org/
* https://mutalyzer.nl/


Interpret variants in other species
----------------------------------

bcfanno designed to annotate VCFs with suitable databases, not restrict to human variants. However the interpret rules may vary from different labs, there is no recommended strategy.



