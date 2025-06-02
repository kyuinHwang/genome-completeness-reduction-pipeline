mkdir ./GffTsv
python ./scripts/join_gff_ko.py --gff=./examples/Gff/GB_GCA_016462095.1.gff --ko=./examples/KO/GB_GCA_016462095.1.ko --output_tsv=./GffTsv/GB_GCA_016462095.1.gff.ko.tsv
python ./scripts/join_gff_ko.py --gff=./examples/Gff/GB_GCA_903959665.1.gff --ko=./examples/KO/GB_GCA_903959665.1.ko --output_tsv=./GffTsv/GB_GCA_903959665.1.gff.ko.tsv
python ./scripts/join_gff_ko.py --gff=./examples/Gff/GB_GCA_947451735.1.gff --ko=./examples/KO/GB_GCA_947451735.1.ko --output_tsv=./GffTsv/GB_GCA_947451735.1.gff.ko.tsv

mkdir ./ReducedKO
python ./scripts/reduce_genome_by_completeness.py --input_tsv=./GffTsv/GB_GCA_016462095.1.gff.ko.tsv --output=./ReducedKO/GB_GCA_016462095.1.reduce.ko --original_completeness=89.8 --target_completeness=48.3
python ./scripts/reduce_genome_by_completeness.py --input_tsv=./GffTsv/GB_GCA_903959665.1.gff.ko.tsv --output=./ReducedKO/GB_GCA_903959665.1.reduce.ko --original_completeness=91.88 --target_completeness=17.2
python ./scripts/reduce_genome_by_completeness.py --input_tsv=./GffTsv/GB_GCA_947451735.1.gff.ko.tsv --output=./ReducedKO/GB_GCA_947451735.1.reduce.ko --original_completeness=54.17 --target_completeness=25.9
