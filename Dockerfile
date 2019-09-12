FROM ubuntu:16.04

# set the data folder as an env var that we can use in scripts
ENV DATA /home/ubuntu/data

# add the local perl libraries to the container and its library search path
ADD ./lib /home/ubuntu/brassica-snps/lib
ENV PERL5LIB "${PERL5LIB}:/home/ubuntu/brassica-snps/lib"

# add the shell and perl script folder to the container and its executable search path
ADD ./script /home/ubuntu/brassica-snps/script
ENV PATH "${PATH}:/home/ubuntu/brassica-snps/script"

# add mart location for ensembl plants
ADD ./conf/biomart /usr/local/src/biomart-perl/conf
ADD ./conf/circos /usr/local/src/circos/conf

# install apt-get packages
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
	sqlite \
	ncbi-blast+ \
	default-jre \
	perl-doc \
	vcftools \
	circos

# install perl packages
RUN cpanm -n \
	Log::Log4perl \
	DBIx::Class \
	DBD::SQLite \
	Bio::Phylo \
	Bio::DB::NCBIHelper \
	Bio::SeqIO \
	Bio::Tools::Run::StandAloneBlastPlus \
	Bio::Tools::Run::RemoteBlast \
	GO::TermFinder \
	Convert::Color \
	SVG \
	git://github.com/naturalis/biomart-perl

# install R packages
RUN R -e "install.packages(c('QTLseqr', 'dplyr', 'ggplot2'), dependencies=T, repos='http://cran.rstudio.com/')"

# all you have to do is this to install GATK because Java is the future
RUN curl -o /usr/local/src/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef'
RUN tar -xvjf /usr/local/src/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
RUN ln -s /usr/local/src/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar /usr/local/src/GenomeAnalysisTK.jar

# and you also have to do this because we need both versions
RUN curl -L -o /usr/local/src/gatk-4.1.3.0.zip https://github.com/broadinstitute/gatk/releases/download/4.1.3.0/gatk-4.1.3.0.zip
RUN unzip /usr/local/src/gatk-4.1.3.0.zip
RUN ln -s /usr/local/src/gatk-4.1.3.0/gatk /usr/local/bin/gatk

# build me as `docker build -t naturalis/brassica-snps .`
# push me as:
# 1. docker login (attempt online through password manager)
# 2. docker push naturalis/brassica-snps

# run me as `docker run -it -v /Users/rutger.vos/deleteme:/home/ubuntu/data naturalis/brassica-snps`
CMD [ "/bin/bash" ]