# Variant types

For genetic variants located in the transcript regions, bcfanno will interpret the amino acid changes and the variant types (**VarType**). One major challange encountered with interpreting VarType is multi transcripts may be encoded by one gene, several VarTypes will be exported and separated with "|" in such a record. Beside VarType, any other transcription-related tags (records), like **HGVSnom**, will also be generated in such a format, which multi transcription-related information seperated with "|" and exported in order.



Supported VarTypes are summaried below.

* [*Intron*](http://www.sequenceontology.org/miso/current_release/term/SO:0001627),  genetic variant in intron;
* [*Noncoding*](http://www.sequenceontology.org/miso/current_release/term/SO:0001792),  genetic variant in the exonic region of a noncoding transcription;
* [*Utr5*](http://www.sequenceontology.org/miso/current_release/term/SO:0001623),  genetic variant in the 5' untranslated region (except splice sites);
* [*Utr3*](http://www.sequenceontology.org/miso/current_release/term/SO:0001624),  genetic variant in the 3' untranslated region (except splice sites);
* [*Synonymous*](http://www.sequenceontology.org/miso/current_release/term/SO:0001819),  genetic variants in coding sequence, resulting in no change to encoded amino acid;
* [*Missense*](http://www.sequenceontology.org/miso/current_release/term/SO:0001583),  genetic variants in coding sequence, resulting in a different amino acid;
* [*Nonsense*](http://www.sequenceontology.org/miso/current_release/term/SO:0001587),  genetic variants in coding sequence, resulting in a premature stop codon, leading to a shortened transcript, also known as *stop-gain*;
* [*InframeInsertion*](http://www.sequenceontology.org/miso/current_release/term/SO:0001821),  insert bases in the coding sequence, resulting in no reading frame shift;
* [*InframeDeletion*](http://www.sequenceontology.org/miso/current_release/term/SO:0001822),  delete bases in the coding sequence, resulting in no reading frame shift;
* [*InframeDelins*](http://www.sequenceontology.org/miso/current_release/term/SO:0001820),  delete and insert bases in the coding sequence, resulting in no reading frame shift;
* [*Frameshift*](http://www.sequenceontology.org/miso/current_release/term/SO:0001589),  delete or insert bases in the coding sequence, resulting in reading frame shift, leading to a incomplete transcript;
* [*StopLost*](http://www.sequenceontology.org/miso/current_release/term/SO:0001578),  genetic variants in the stop codon, resulting in an elongated transcript;
* [*StopRetained*](http://www.sequenceontology.org/miso/current_release/term/SO:0001567),  genetic variants in the stop codon but terminator remains;
* [*SpliceSite*](http://www.sequenceontology.org/miso/current_release/term/SO:0001630),  genetic variants in the region of splice site, either within 1-3 bases of the exon or 3-8 bases of the intron, please notice that 1-2 bases of the intron will be interpret as *SpliceDonor* or *SpliceAcceptor*; 
* [*SpliceDonor*](http://www.sequenceontology.org/miso/current_release/term/SO:0001575),  genetic variants in the 2 base region at the 5' end of an intron;
* [*SpliceAcceptor*](http://www.sequenceontology.org/miso/current_release/term/SO:0001574), genetic variants in the 2 base region at the 3' end of an intron;
* *Complex*,  large variants overlapped coding and noncoding regions, or influence more than one gene;
* *NoCall*,  genetic variants in coding region, and the alternative allele is the same with the base in the transcript sequence, the *NoCall* type come from the inconsistance between genome reference and transcription reference, and usually come with high allele frequency in population database;
* *Unknown*, genetic variants in coding region but could not be interpreted by bcfanno, usually account for program bugs.

For any other variants not annotated with bcfanno (empty or no VarType tag) could be interpret as intergenic type.



# HGVS nomenclature













