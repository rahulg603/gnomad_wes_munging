import argparse
import hail as hl
import dxpy
import pyspark
import re
from gnomad.utils import vep


MNT_DIR = 'file:///mnt/project/Bulk/'
MT_META_GNOMAD = 'ukb24068_c21_b0_v1.vcf.gz'
MT_META_UKB = 'ukb23148_c21_b0_v1.vcf.gz'
FEATURE_FIELD_TO_TYPE = {"Transcript": ("transcript_consequences", "transcript_feature_id"),
                         "RegulatoryFeature": ("regulatory_feature_consequences", "regulatory_feature_id"),
                         "MotifFeature": ("motif_feature_consequences", "motif_feature_id"),
                         "Intergenic": ("intergenic_consequences", "")}
CHROM = [str(x) for x in range(1, 23)] + ['X', 'Y']
LOF_CRITERIA = lambda x: (x.lof == 'HC') & (hl.len(x.lof_flags) == 0)

# see Loftee docs here for information on where these came from: 
# https://github.com/broadinstitute/gnomad-browser/blob/main/browser/src/VariantPage/Loftee.tsx
PROBLEM_PREFIXES = ['PHYLOCSF_', 'NON_CAN_SPLICE', 'NAGNAG_SITE', 'SINGLE_EXON',
                    'GERP_DIST',
                    'BP_DIST',
                    'DIST_FROM_LAST_EXON',
                    '50_BP_RULE',
                    'PERCENTILE:',
                    'ANN_ORF:',
                    'MAX_ORF:',
                    'INTRON_SIZE:',
                    'DE_NOVO_DONOR']
# flags and info -- look into LOFTEE docs
LOFTEE_FILTERS = ['END_TRUNC', 'INCOMPLETE_CDS', 
                  'EXON_INTRON_UNDEF', 'SMALL_INTRON',
                  'ANC_ALLELE', 'NON_DONOR_DISRUPTING',
                  'NON_ACCEPTOR_DISRUPTING', 'RESCUE_DONOR',
                  'RESCUE_ACCEPTOR', 'GC_TO_GT_DONOR', 
                  '5UTR_SPLICE', '3UTR_SPLICE']


def parse_positions(expr):
    raise NotImplementedError('Position parsing is not yet implemented.')
    return 1, 1


def parse_score(expr):
    exprnew = expr.replace('\\)$', '').split('\\(', 2)
    expr_pred = hl.or_missing((expr.length() > 0) & (hl.len(exprnew) == 2), exprnew[0])
    expr_score = hl.or_missing((expr.length() > 0) & (hl.len(exprnew) == 2) & (hl.is_defined(exprnew[1])), hl.float64(exprnew[1]))
    return expr_pred, expr_score


def format_vep(mt, mt_meta, temp_dir, skip_vep_check, fix_fields):
    """
    We parse tradiitonal VEP outputs based on the VCF metadata.
    We also note that in some instances a few new fields are present.
    """
    str_vep_format = mt_meta['info']['vep']['Description']
    format_info = re.search('.+Format: (.+)$', str_vep_format)[1].lower()
    format_info_list = format_info.split('|')
    
    # there is an edge case where some comma-delimited strings are inserted, causing a parsing issue.
    # we solve this by remerging vep, replacing the commas with semicolons, and then respltting twice
    if fix_fields:
        mt = mt.annotate_rows(vep_merged = hl.literal(',').join(mt.info.vep))
        mt = mt.annotate_rows(vep_merged = mt.vep_merged.replace(',(?=' + '|'.join(PROBLEM_PREFIXES) + ')', ';'))
        mt = mt.annotate_rows(vep_merged = mt.vep_merged.replace('(?<=' + '|'.join(LOFTEE_FILTERS) + ')' + ',(?=' + '|'.join(LOFTEE_FILTERS) + ')', ';'))
        mt = mt.annotate_rows(vep_split = mt.vep_merged.split(',').map(lambda x: x.split('\\|')))
        mt = mt.drop('vep_merged')
    else:
        mt = mt.annotate_rows(vep_split = mt.info.vep.map(lambda x: x.split('\\|')))

    # previously we removed records with only 1 element in VEP after splitting;
    # we now require that the above approach solves all of the issues
    # mt = mt.annotate_rows(vep_split = hl.filter(lambda y: hl.len(y) != 1, mt.vep_split))
    if not skip_vep_check:
        mt = mt.annotate_rows(vep_split_ct_non45 = hl.sum(mt.vep_split.map(lambda x: hl.if_else(hl.len(x) != 45, 1, 0))))
        mt_filtered_wrong_len = mt.filter_rows(mt.vep_split_ct_non45 != 0)
        mt_filtered_wrong_len = mt_filtered_wrong_len.repartition(200).checkpoint(f'{temp_dir}checking_incorrect_vep_lengths.mt', overwrite=True)
        count_wrong_vep = mt_filtered_wrong_len.count_rows()
        if count_wrong_vep > 0:
            print(mt_filtered_wrong_len.info.vep.show())
            raise ValueError(f'ERROR: {str(count_wrong_vep)} records found to have incorrect length VEP string.')
        mt = mt.drop('vep_split_ct_non45')
    

    def return_csq_to_struct(veparray):
        """ Reverses the gnomAD VEP tool `vep_struct_to_csq`
        Does not currently deal with domains as this is not currently needed.
        Does not currently do float conversion.
        Does not currently deal with positions or scores.
        """
        vep_struct = hl.struct()
        for_int32 = ['allele_num', 'distance', 'gene_pheno', 'hgvs_offset', 'minimised', 'strand', 'tsl', 'motif_pos']
        #for_float64 = ['motif_score_change']
        # Global modifications
        veparray = veparray.map(lambda element: element.annotate(**{'variant_allele': element.allele,
                                                                    'consequence_terms': element.consequence.split('&')}
                                                                 ).drop('allele', 'consequence'))
        # Type changes
        veparray = veparray.map(lambda element: element.annotate(**{x: hl.or_missing((element[x].length() > 0) & (hl.is_defined(element[x])), hl.int32(element[x])) for x in for_int32}))
        #veparray = veparray.map(lambda element: element.annotate(**{x: hl.or_missing((element[x].length() > 0) & (hl.is_defined(element[x])), hl.float64(element[x])) for x in for_float64}))
        # Per feature modifications
        for feature, (feature_field, feature_id) in FEATURE_FIELD_TO_TYPE.items():
            this_item = hl.filter(lambda x: x.feature_type == feature, veparray)
            if feature != 'Intergenic':
                this_item = this_item.map(lambda element: element.annotate(**{feature_id: element.feature}).drop('feature'))
            if feature == 'Transcript':
                this_item = this_item.map(lambda element: element.annotate(**{'canonical': hl.if_else(element.canonical == 'YES', 1, 0),
                                                                              'polyphen_prediction': parse_score(element.polyphen)[0],
                                                                              #'polyphen_score': parse_score(element.polyphen)[1],
                                                                              'sift_prediction': parse_score(element.sift)[0]}
                                                                              #'sift_score': parse_score(element.sift)[1]}
                                                                          ).rename({'ensp':'protein_id',
                                                                                    'gene':'gene_id',
                                                                                    'symbol':'gene_symbol',
                                                                                    'symbol_source':'gene_symbol_source'}))
            else: # dropping processed fields from other feature types
                this_item = this_item.map(lambda element: element.drop('polyphen','sift','protein_position','cds_position','cdna_position'))
            vep_struct = vep_struct.annotate(**{feature_field: this_item})
        
        return vep_struct


    csq = mt.vep_split.map(lambda x: hl.struct(**{y: x[idy] for idy, y in enumerate(format_info_list)}))
    mt = mt.annotate_rows(vep = return_csq_to_struct(csq))
    mt = mt.drop('vep_split')
    mt = mt.annotate_rows(info = mt.info.drop('vep'))
    return mt


def generate_lof_mt(args, data_path, mid_str, batch_pref, checkpoint_mt_path, temp_dir, apply_filters_after):
    mt_meta = hl.get_vcf_metadata(f'{data_path}/{MT_META_UKB if args.use_default_exomes else MT_META_GNOMAD}')
    if args.fix_fields_using_internal_method:
        mt = hl.import_vcf(f'{data_path}/*{mid_str}_b{batch_pref}*.vcf.gz', 
                           reference_genome='GRCh38', 
                           force_bgz=(not args.use_default_exomes))
    else:
        mt = hl.import_vcf(f'{data_path}/*{mid_str}_b{batch_pref}*.vcf.gz', 
                           reference_genome='GRCh38', 
                           force_bgz=(not args.use_default_exomes),
                           find_replace=(',(?=' + '|'.join(PROBLEM_PREFIXES + LOFTEE_FILTERS) + ')',';'))
    
    mt_formatted = format_vep(mt, mt_meta, temp_dir=temp_dir, skip_vep_check=args.skip_vep_check, fix_fields=args.fix_fields_using_internal_method)
    if not apply_filters_after:
        mt_formatted = mt_formatted.filter_rows(hl.len(mt_formatted.filters) == 0)
    mt_csq = vep.process_consequences(mt_formatted)

    # Now filter to LOF variants
    mt_csq = mt_csq.annotate_rows(worst_csq_gene = mt_csq.vep.worst_csq_by_gene_canonical
                  ).explode_rows('worst_csq_gene')
    if not apply_filters_after:
        mt_csq_lof = mt_csq.filter_rows(LOF_CRITERIA(mt_csq.worst_csq_gene))
    else:
        mt_csq_lof = mt_csq.filter_rows(hl.literal(vep.LOF_CSQ_SET).contains(mt_csq.worst_csq_gene.most_severe_consequence))
    
    if args.allow_all_ids:
        mt_csq_lof = mt_csq_lof.annotate_rows(variant_site = hl.or_missing(mt_csq_lof.worst_csq_gene.hgvsc.length() > 0, 
                                                                           mt_csq_lof.worst_csq_gene.hgvsc.split(':', 2)[1]))
        mt_csq_lof = mt_csq_lof.key_rows_by('locus', 'alleles', 'variant_site').distinct_by_row()
    else:
        mt_csq_lof = mt_csq_lof.filter_rows(mt_csq_lof.worst_csq_gene.gene_symbol_source == 'HGNC')

    # Checkpoint the table
    mt_csq_lof = mt_csq_lof.repartition(200).checkpoint(checkpoint_mt_path, overwrite=True)
    return mt_csq_lof


def generate_gene_lof_summary(mt):
    """ Adapted from https://github.com/broadinstitute/gnomad_lof/blob/e09d9fdc13d43c8ac996c93735ed206ba41ed473/constraint/gene_lof_matrix.py#L83
    """
    ht = mt.annotate_rows(no_lofs=hl.agg.count_where((mt.defined_sites > 0) & (mt.num_hom_alt + mt.num_het == 0)),
                          obs_het_lof=hl.agg.count_where((mt.num_het > 0) & (mt.num_hom_alt == 0)),
                          obs_hom_lof=hl.agg.count_where(mt.num_hom_alt > 0),
                          defined=hl.agg.count_where(mt.defined_sites > 0),
                          ).rows()
    ht = ht.annotate(p=1 - hl.sqrt(hl.float64(ht.no_lofs) / ht.defined))
    ht = ht.annotate(exp_hom_lof=ht.defined * ht.p * ht.p)
    return ht.annotate(oe=ht.obs_hom_lof / ht.exp_hom_lof)


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def chunk_vector(vec, chunk_size):
    thislen = len(vec)
    if thislen <= chunk_size:
        return [vec], [str(0)]
    else:
        vec_out = list(chunks(vec, chunk_size))
        return vec_out, [str(x) for x in range(0, len(vec_out))]


parser = argparse.ArgumentParser()
parser.add_argument('--use-default-exomes', help='If enabled, will import UKB processed 450k exomes. If not enabled, will import gnomAD exomes.',
                    action='store_true')
parser.add_argument('--specific-chr', type=str, default='', help='Only import a specific chromosome.')
parser.add_argument('--dx-init', help='DNANexus SQL database.')
parser.add_argument('--temp-dir', help='Relative path to temporary file directory.')
parser.add_argument('--export-prefix', help='Relative path to exports. Will output the per-sample tsv and a checkpointed mt.')
parser.add_argument('--batch-prefix', type=str, help='Troubleshooting prefix for testing purposes. Limits the number of batches imported.')
parser.add_argument('--allow-all-ids', help='Allows all IDs. Will not try to select a particular ID type. If disabled, filters to EntrezGene.')
parser.add_argument('--overwrite-filtered-table', help='If true, will overwrite the filtered table.', action='store_true')
parser.add_argument('--overwrite-after-filtering', action='store_true', help='If enabled, will overwrite after filtering.')
parser.add_argument('--skip-vep-check', help='If enabled, will avoid the VEP length check. Adds speed but will lead to potentially strange parsing issues during export.', action='store_true')
parser.add_argument('--fix-fields-using-internal-method', help='If enabled, will split using regex. Default is to fix fields on import.', action='store_true')
parser.add_argument('--assume-all-files-exist', action='store_true', help='If enabled, will assume MTs all exist.')
parser.add_argument('--apply-filters-post-checkpoint', action='store_true', help='If enabled, will apply filters after checkpoint. Two filters are applied -- removal of anything with FT fields and filtering to LOFTEE interpreted LOF alleles. Enable this as a band-aid solution to avoid exporting all the exomes again but apply filters.')
parser.add_argument('--avoid-generating-full-lof-table', action='store_true', help='If enabled, will skip generation of full LoF table. This amounts to a union of all credible LoFs.')
parser.add_argument('--merging-chunk-size', type=int, default=3, help='Will use these chunk sizes to speed up merging.')


if __name__ == '__main__':
    args = parser.parse_args()
    my_database = dxpy.find_one_data_object(name=args.dx_init.lower())['id']
    sc = pyspark.SparkContext()
    spark = pyspark.sql.SparkSession(sc)
    temp_dir = f'dnax://{my_database}/{args.temp_dir}/'
    hl.init(sc=sc, tmp_dir=temp_dir)
    suffpath = 'Exome sequences/Population level exome OQFE variants, pVCF format - interim 450k release' if args.use_default_exomes else 'Exome sequences_Alternative exome processing/Exome variant call files (gnomAD) (VCFs)'
    data_path = MNT_DIR + suffpath
    tf_overwrite_next = args.overwrite_filtered_table or args.overwrite_after_filtering

    CHROM_VEC = [args.specific_chr] if len(args.specific_chr) != 0 else CHROM
    mts_csq_lof = []
    mts_per_gene = []
    s_list = []
    batch_pref = '' if args.batch_prefix is None else args.batch_prefix
    batch_pref_this = '' if batch_pref == '' else f'_b{batch_pref}'
    for chr in CHROM_VEC:
        mid_str = f'c{chr}'
        print(f'Loading {mid_str}...')

        checkpoint_mt_path = f'dnax://{my_database}/{args.export_prefix}_worstcsqlof_filtered_{mid_str}{batch_pref_this}.mt'
        if not args.assume_all_files_exist:
            if (not args.overwrite_filtered_table) & hl.hadoop_is_file(f'{checkpoint_mt_path}/_SUCCESS'):
                mt_csq_lof = hl.read_matrix_table(checkpoint_mt_path)
                print(f'Loaded {mid_str} LoF table from disk.')
            else:
                mt_csq_lof = generate_lof_mt(args, data_path, mid_str, 
                                             batch_pref, checkpoint_mt_path, temp_dir,
                                             args.apply_filters_post_checkpoint)

            if args.apply_filters_post_checkpoint:
                print(f'Assuming {mid_str} LoF table needs further filtering post-checkpoint...')
                mt_csq_lof = mt_csq_lof.filter_rows(hl.len(mt_csq_lof.filters) == 0).repartition(200)
                mt_csq_lof = mt_csq_lof.filter_rows(LOF_CRITERIA(mt_csq_lof.worst_csq_gene))
            mts_csq_lof.append(mt_csq_lof)

        # Now, per gene, compute the number of loss of function variants found per individual
        per_gene_mt_path = f'dnax://{my_database}/{args.export_prefix}_worstcsqlof_filtered_pergene_{mid_str}{batch_pref_this}.mt'
        if args.assume_all_files_exist:
            per_gene_group = hl.read_matrix_table(per_gene_mt_path)
        else:
            if (not tf_overwrite_next) & hl.hadoop_is_file(f'{per_gene_mt_path}/_SUCCESS'):
                print(f'Loaded {mid_str} per-gene LoF counts table from disk.')
                per_gene_group = hl.read_matrix_table(per_gene_mt_path)
            else:
                per_gene_group = mt_csq_lof.group_rows_by(mt_csq_lof.worst_csq_gene.gene_symbol
                                          ).aggregate(num_hom_alt = hl.agg.count_where(mt_csq_lof.GT.is_hom_var()),
                                                      num_het = hl.agg.count_where(mt_csq_lof.GT.is_het()),
                                                      defined_sites=hl.agg.count_where(hl.is_defined(mt_csq_lof.GT)))
                per_gene_group = per_gene_group.repartition(10).checkpoint(f'dnax://{my_database}/{args.export_prefix}_worstcsqlof_filtered_pergene_{mid_str}{batch_pref_this}.mt', overwrite=tf_overwrite_next)
        mts_per_gene.append(per_gene_group)
        s_list.append(per_gene_group.s.collect())

        if not args.assume_all_files_exist:
            # Produce summary data
            flat_file_summ_path = f'dnax://{my_database}/{args.export_prefix}_worstcsqlof_filtered_pergene_summarystats_{mid_str}{batch_pref_this}.tsv.bgz'
            if tf_overwrite_next or (not hl.hadoop_is_file(flat_file_summ_path)):
                print(f'Generating {mid_str} per-gene LoF summaries...')
                ht = generate_gene_lof_summary(per_gene_group)
                ht.export(flat_file_summ_path)
            
            # Subset to 0, 1, 2 for a given gene
            flat_file_tsv_path = f'dnax://{my_database}/{args.export_prefix}_worstcsqlof_filtered_pergene_gt_{mid_str}{batch_pref_this}.tsv.bgz'
            if tf_overwrite_next or (not hl.hadoop_is_file(flat_file_tsv_path)):
                print(f'Generating {mid_str} per-gene LoF counts...')
                per_gene_lofgt = per_gene_group.select_entries(GT_LOF = hl.case(
                                                                        ).when(per_gene_group.num_hom_alt > 0, 2
                                                                        ).when(per_gene_group.num_het > 0, 1
                                                                        ).default(0))
                per_gene_lofgt.GT_LOF.export(flat_file_tsv_path)

    if len(args.specific_chr) == 0:
        print("Per chromosome sample sizes (should be equal):")
        print({'chr'+str(idx):len(x) for idx, x in enumerate(s_list)})
        vec_shared = s_list[0]
        if len(CHROM_VEC) > 1:
            for idx, chr in enumerate(CHROM_VEC[1:]):
                vec_shared = [vec_shared[x] if vec_shared[x] == s_list[idx+1][x] else 'RG_INCORRECT_FOR_RM' for x in range(0, len(vec_shared))]
        vec_shared_fin = [x for x in vec_shared if x != 'RG_INCORRECT_FOR_RM']
        print(f'Filtered sample vector from {str(len(vec_shared))} to {str(len(vec_shared_fin))}.')

        final_mt_path = f'dnax://{my_database}/{args.export_prefix}_worstcsqlof_filtered_allchr{batch_pref_this}.mt'
        
        if not args.avoid_generating_full_lof_table:
            mts_csq_lof_f = []
            for idx, chr in enumerate(CHROM_VEC):
                this_mt = mts_csq_lof[idx]
                mts_csq_lof_f.append(this_mt.filter_cols(hl.literal(vec_shared_fin).contains(this_mt.s)))
            full_mt_csq_lof = hl.MatrixTable.union_rows(*mts_csq_lof_f).repartition()
            full_mt_csq_lof.write(final_mt_path, overwrite=True)
        
        mts_per_gene_f = []
        for idx, chr in enumerate(CHROM_VEC):
            this_mt = mts_per_gene[idx]
            mts_per_gene_f.append(this_mt.filter_cols(hl.literal(vec_shared_fin).contains(this_mt.s)))
        
        mts_per_gene_intermediate = []
        chunked_mts, chunked_idx = chunk_vector(mts_per_gene_f, args.merging_chunk_size)
        for mtchunk, idx in zip(chunked_mts, chunked_idx):
            this_temp_path = f'{temp_dir}temp_chunked_per_gene_{idx}.mt'
            this_mrg_per_gene = hl.MatrixTable.union_rows(*mtchunk
                                             ).checkpoint(this_temp_path, overwrite=True)
            mts_per_gene_intermediate.append(this_mrg_per_gene)
        
        full_mt_per_gene_group = hl.MatrixTable.union_rows(*mts_per_gene_intermediate
                                              ).checkpoint(f'dnax://{my_database}/{args.export_prefix}_worstcsqlof_filtered_pergene_allchr{batch_pref_this}.mt', overwrite=True)

        full_ht = generate_gene_lof_summary(full_mt_per_gene_group)
        full_ht.export(f'dnax://{my_database}/{args.export_prefix}_worstcsqlof_filtered_pergene_summarystats_allchr{batch_pref_this}.tsv.bgz')

        full_per_gene_lofgt = full_mt_per_gene_group.select_entries(GT_LOF = hl.case(
                                                                              ).when(full_mt_per_gene_group.num_hom_alt > 0, 2
                                                                              ).when(full_mt_per_gene_group.num_het > 0, 1
                                                                              ).default(0))
        full_per_gene_lofgt.GT_LOF.export(f'dnax://{my_database}/{args.export_prefix}_worstcsqlof_filtered_pergene_gt_allchr{batch_pref_this}.tsv.bgz')
