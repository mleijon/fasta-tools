#!/usr/bin/python


class SeqFeatures:
    def __init__(self, seq):
        self.seq_locus = seq.find('INSDSeq_locus').text
        self.seq_length = seq.find('INSDSeq_length').text
        self.seq_strandedness = seq.find('INSDSeq_strandedness').text
        self.seq_moltype = seq.find('INSDSeq_moltype').text
        self.seq_topology = seq.find('INSDSeq_topology').text
        self.seq_division = seq.find('INSDSeq_division').text
        self.seq_updateDate = seq.find('INSDSeq_update-date').text
        self.seq_createDate = seq.find('INSDSeq_create-date').text
        self.seq_definition = seq.find('INSDSeq_definition').text
        self.seq_primaryAccession = seq.find('INSDSeq_primary-accession').text
        self.seq_accessionVersion = seq.find('INSDSeq_accession-version').text
        self.seq_OtherSeqIds = list()
        for item in seq.findall('INSDSeq_other-seqids/INSDSeqid'):
            self.seq_OtherSeqIds.append(item.text)
        self.seq_keywords = list()
        if not seq.find('INSDSeq_project') is None:
            self.seq_project = seq.find('INSDSeq_project').text
        for item in seq.findall('INSDSeq_keywords/INSDKeyword'):
            self.seq_keywords.append(item.text)
        self.seq_source = seq.find('INSDSeq_source').text
        self.seq_organism = seq.find('INSDSeq_organism').text
        self.seq_taxonomy = seq.find('INSDSeq_taxonomy').text
        self.seq_references = list()
        reference = dict()
        for item in seq.findall('INSDSeq_references/INSDReference'):
            reference['reference'] = item.find('INSDReference_reference').text
            reference['reference_position'] = \
                item.find('INSDReference_position').text
            if not item.find('INSDReference_consortium') is None:
                reference['reference_consortium'] = \
                    item.find('INSDReference_consortium').text
            reference['reference_title'] = item.find('INSDReference_title').text
            reference['reference_journal'] = \
                item.find('INSDReference_journal').text
            reference['reference_title'] = item.find('INSDReference_title').text
            if not item.find('INSDReference_authors') is None:
                authors = list()
                for author in item.findall('INSDReference_authors/INSDAuthor'):
                    authors.append(author.text)
                reference['reference_authors'] = authors
            self.seq_references.append(reference)
            reference = dict()
        self.seq_comment = seq.find('INSDSeq_comment').text
        self.seq_sequence = seq.find('INSDSeq_sequence').text


class FivepUTR:
    def __init__(self, ft):
        self.feature_key = ft.find('INSDFeature_key').text
        self.feature_location = ft.find('INSDFeature_location').text
        interval = dict()
        intervallst = list()
        for item in ft.findall('INSDFeature_intervals/INSDInterval'):
            if not item.find('INSDInterval_from') is None:
                interval['interval_from'] = item.find('INSDInterval_from').text
                interval['interval_to'] = item.find('INSDInterval_to').text
                interval['interval_acc'] = item.find(
                    'INSDInterval_accession').text
            elif not item.find('INSDInterval_point') is None:
                interval['interval_from'] = item.find('INSDInterval_point').text
                interval['interval_to'] = item.find('INSDInterval_point').text
                interval['interval_acc'] = item.find(
                    'INSDInterval_accession').text
            else:
                interval['interval_from'] = ""
                interval['interval_to'] = ""
                interval['interval_acc'] = ""
            intervallst.append(interval)
            interval = dict()
        self.feature_intervals = intervallst


class Source(FivepUTR):
    def __init__(self, ft):
        super().__init__(ft)
        qualifiers = dict()
        for qualifier in ft.findall('INSDFeature_quals/INSDQualifier'):
            if not qualifier.find('INSDQualifier_name') is None:
                name = qualifier.find('INSDQualifier_name').text
            else:
                continue
            if not qualifier.find('INSDQualifier_value') is None:
                value = qualifier.find('INSDQualifier_value').text
            else:
                value = ""
            qualifiers[name] = value
        self.qualifiers = qualifiers


class SeqGene(Source):
    def __init__(self, ft):
        super().__init__(ft)


class CDS(Source):
    def __init__(self, ft):
        super().__init__(ft)
        if not ft.find('INSDFeature_operator') is None:
            self.feature_operator = ft.find('INSDFeature_operator').text


class MatPeptide(CDS):
    def __init__(self, ft):
        super().__init__(ft)


class ThreepUTR(FivepUTR):
    def __init__(self, ft):
        super().__init__(ft)


class MiscFeature(Source):
    def __init__(self, ft):
        super().__init__(ft)


class MiscRNA(Source):
    def __init__(self, ft):
        super().__init__(ft)


class StemLoop(FivepUTR):
    def __init__(self, ft):
        super().__init__(ft)


class Regulatory(Source):
    def __init__(self, ft):
        super().__init__(ft)


def extract_seq(xmlroot):
    outfile = open(ARGS.o, 'w')
    for sequence in xmlroot.iter('INSDSeq'):
        gene_name = ""
        cds_name = ""
        seq_name = ""
        mp_name = ""
        seq = ""
        int_from = set()
        int_to = set()
        seqfeatures = SeqFeatures(sequence)
        for feature in sequence.findall(
                'INSDSeq_feature-table/INSDFeature'):
            acc = seqfeatures.seq_primaryAccession
            organism = seqfeatures.seq_organism.replace(' ', '_')
            if ARGS.t == 's':
                if feature.find('INSDFeature_key').text == 'source':
                    source = Source(feature)
                    for interval in source.feature_intervals:
                        int_from.add(int(interval['interval_from']))
                        int_to.add(int(interval['interval_to']))
                    min_int = min(int_from)
                    max_int = max(int_to)
                    seq = seqfeatures.seq_sequence[min_int - 1: max_int]
                seq_name = '>' + acc + ';' + organism + '\n'
            if ARGS.t == 'g':
                if feature.find('INSDFeature_key').text == 'gene':
                    gene = SeqGene(feature)
                    if 'gene' not in gene.qualifiers.keys():
                        continue
                    elif gene.qualifiers['gene'].casefold() not in \
                            feature_names:
                        continue
                    gene_name = gene.qualifiers['gene']
                    for interval in gene.feature_intervals:
                        int_from.add(int(interval['interval_from']))
                        int_to.add(int(interval['interval_to']))
                    min_int = min(int_from)
                    max_int = max(int_to)
                    seq = seqfeatures.seq_sequence[min_int - 1: max_int]
                seq_name = '>' + acc + ';' + organism + ';gene=' + gene_name\
                           + '\n'
            if ARGS.t == 'c':
                if feature.find('INSDFeature_key').text == 'CDS':
                    cds = CDS(feature)
                    if 'product' not in cds.qualifiers.keys():
                        continue
                    elif cds.qualifiers['product'].casefold() not in \
                            feature_names:
                        continue
                    cds_name = cds.qualifiers['product']
                    if ARGS.s == 'a':
                        seq = cds.qualifiers['translation']
                    else:
                        for interval in cds.feature_intervals:
                            int_from.add(int(interval['interval_from']))
                            int_to.add(int(interval['interval_to']))
                        min_int = min(int_from)
                        max_int = max(int_to)
                        seq = seqfeatures.seq_sequence[min_int - 1: max_int]
                seq_name = '>' + acc + ';' + organism + ';CDS=' + cds_name \
                           + '\n'
            if ARGS.t == 'm':
                if feature.find('INSDFeature_key').text == 'mat_peptide':
                    mat_peptide = MatPeptide(feature)
                    if 'product' not in mat_peptide.qualifiers.keys():
                        continue
                    elif mat_peptide.qualifiers['product'].casefold() not in \
                            feature_names:
                        continue
                    mp_name = mat_peptide.qualifiers['product']
                    if ARGS.s == 'a':
                        seq = mat_peptide.qualifiers['peptide']
                    else:
                        for interval in mat_peptide.feature_intervals:
                            int_from.add(int(interval['interval_from']))
                            int_to.add(int(interval['interval_to']))
                        min_int = min(int_from)
                        max_int = max(int_to)
                        seq = seqfeatures.seq_sequence[min_int - 1: max_int]
                seq_name = '>' + acc + ';' + organism + ';peptide=' + mp_name \
                           + '\n'
            if ARGS.t == '5':
                if feature.find('INSDFeature_key').text == '5\'UTR':
                    fivep_utr = FivepUTR(feature)
                    for interval in fivep_utr.feature_intervals:
                        int_from.add(int(interval['interval_from']))
                        int_to.add(int(interval['interval_to']))
                    min_int = min(int_from)
                    max_int = max(int_to)
                    seq = seqfeatures.seq_sequence[min_int - 1: max_int]
                seq_name = '>' + acc + ';' + organism + ';5\'UTR\n'
            if ARGS.t == '3':
                if feature.find('INSDFeature_key').text == '3\'UTR':
                    threep_utr = ThreepUTR(feature)
                    print(threep_utr.feature_intervals)
                    for interval in threep_utr.feature_intervals:
                        int_from.add(int(interval['interval_from']))
                        int_to.add(int(interval['interval_to']))
                    min_int = min(int_from)
                    max_int = max(int_to)
                    seq = seqfeatures.seq_sequence[min_int - 1: max_int]
                seq_name = '>' + acc + ';' + organism + ';3\'UTR\n'
        if seq != "":
            outfile.write(seq_name + seq + '\n')
    outfile.close()


def extract_info(xmlroot, inff):
    genes = dict()
    cdss = dict()
    mps = dict()
    for sequence in xmlroot.iter('INSDSeq'):
        for feature in sequence.findall(
                'INSDSeq_feature-table/INSDFeature'):
            if feature.find('INSDFeature_key').text == 'gene':
                gene = SeqGene(feature)
                if 'gene' not in gene.qualifiers.keys():
                    continue
                if gene.qualifiers['gene'].casefold() in genes.keys():
                    genes[gene.qualifiers['gene'].casefold()] += 1
                else:
                    genes[gene.qualifiers['gene'].casefold()] = 1
            elif feature.find('INSDFeature_key').text == 'CDS':
                cds = CDS(feature)
                if cds.qualifiers['product'].casefold() in cdss.keys():
                    cdss[cds.qualifiers['product'].casefold()] += 1
                else:
                    cdss[cds.qualifiers['product'].casefold()] = 1
            elif feature.find('INSDFeature_key').text == 'mat_peptide':
                mp = MatPeptide(feature)
                if mp.qualifiers['product'].casefold() in mps.keys():
                    mps[mp.qualifiers['product'].casefold()] += 1
                else:
                    mps[mp.qualifiers['product'].casefold()] = 1
    inff.write('GENES:\n{}\n\nCDS:\n{}\n\nMATURE PEPTIDES:\n{}'.format(
        dict(sorted(genes.items(), key=lambda item: item[1], reverse=True)),
        dict(sorted(cdss.items(), key=lambda item: item[1], reverse=True)),
        dict(sorted(mps.items(), key=lambda item: item[1], reverse=True))))


if __name__ == "__main__":
    import argparse
    import xml.etree.ElementTree as Et

    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-f', type=str, help='input xml-file', required=True)
    PARSER.add_argument('-o', type=str, help='output file', required=True)
    PARSER.add_argument('-t', choices=['s', 'g', 'c', 'm', '5', '3'], type=str,
                        help='type of feature (source[s] = full nt-sequence;'
                             'gene[g]; CDS[c]; mature peptide[m]; 5\'UTR[5];'
                             '3\'UTR[3])', default='s')
    PARSER.add_argument('-n', type=str, help='Feature name', nargs='+', )
    PARSER.add_argument('-s', choices=['a', 'n'], type=str,
                        help='Sequence type', default='n')
    PARSER.add_argument('-l', action="store_true", default=False,
                        help='Switch for output og info-file')
    ARGS = PARSER.parse_args()
    if ARGS.t in ['g', 'c', 'm'] and not ARGS.n:
        exit('No feature name. Use the -n option to give the name of the '
             'selected feature. Note: this can be a list of names separated by'
             'spaces')
    if ARGS.t in ['g', '5', '3'] and ARGS.s == 'a':
        exit('Genes and UTRs can only have nt-sequence type. Leave out the -s'
             ' option or set -s n')
    if ARGS.n:
        feature_names = [item.casefold() for item in ARGS.n]
    tree = Et.parse(ARGS.f)
    root = tree.getroot()
    if ARGS.l:
        info_file = open(ARGS.f.split('.')[0] + '.info', 'w')
        extract_info(root, info_file)
    extract_seq(root)
