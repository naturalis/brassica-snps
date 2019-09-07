FROM ubuntu:16.04
ADD ./lib /brassica-snps/lib
ADD ./script /brassica-snps/script
ADD ./data /brassica-snps/data
RUN apt install cpanminus bwa samtools muscle sqlite r-base graphviz libgd-perl
RUN cpanm -n Log::Log4perl DBIx::Class DBD::SQLite Bio::Phylo Bio::DB::NCBIHelper Bio::SeqIO Bio::Tools::Run::StandAloneBlastPlus GO::TermFinder git://github.com/naturalis/biomart-perl
RUN Rscript -e 'install.packages("QTLseqr", "dplyr", "ggplot2")'
ENV PERL5LIB "${PERL5LIB}:/brassica-snps/lib"
ENV PATH "${PATH}:/brassica-snps/script"
CMD [ "python", "/imgpheno/examples/sticky-traps.py" ]

