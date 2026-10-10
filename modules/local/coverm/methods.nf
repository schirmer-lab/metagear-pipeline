// CoverM names each column "<bam stem> <suffix>"; the modules strip the suffix so the
// published header is the sample. Suffixes as printed by coverm 0.7.0.
def covermColumnSuffix( method ) {
    def suffixes = [
        count:            'Read Count',
        trimmed_mean:     'Trimmed Mean',
        rpkm:             'RPKM',
        tpm:              'TPM',
        covered_bases:    'Covered Bases',
        covered_fraction: 'Covered Fraction',
        length:           'Length',
    ]

    if ( !suffixes.containsKey(method) ) {
        error "coverm method '${method}' has no column suffix in modules/local/coverm/methods.nf; " +
              "run `coverm contig --methods ${method}` and add the header text it prints"
    }

    return suffixes[method]
}
