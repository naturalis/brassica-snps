FROM ubuntu:16.04
ADD ./lib /home/ubuntu/brassica-snps/lib
ADD ./script /home/ubuntu/brassica-snps/script
RUN apt-get update && apt-get install -y \
	bwa \
	cpanminus \
	git-core \
	graphviz \
	mercurial \
	libgd-perl \
	muscle \
	r-base \
	samtools \
	sqlite  
RUN cpanm -n \
	Log::Log4perl \
	DBIx::Class \
	DBD::SQLite \
	Bio::Phylo \
	Bio::DB::NCBIHelper \
	Bio::SeqIO \
	Bio::Tools::Run::StandAloneBlastPlus \
	GO::TermFinder \
	git://github.com/naturalis/biomart-perl
RUN R -e "install.packages(c('QTLseqr', 'dplyr', 'ggplot2'), dependencies=T, repos='http://cran.rstudio.com/')"
ENV PERL5LIB "${PERL5LIB}:/home/ubuntu/brassica-snps/lib"
ENV PATH "${PATH}:/home/ubuntu/brassica-snps/script"
CMD [ "python", "/imgpheno/examples/sticky-traps.py" ]

