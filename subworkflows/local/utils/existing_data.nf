/* --- Utility Functions --- */

/**
 * Handle the common pattern of creating existing files channel and splitting input
 * into items that need processing vs items that already exist
 */
def createExistingDirChannel(param_path, file_pattern, file_suffix, transform_fn = null) {

    // Create existing data channel
    if ( !param_path ) {
        return Channel.empty()
    } else {
        def files_path = file(param_path)

        if ( !files_path.isDirectory() ) {
            error "The provided path is not a directory: ${files_path}"
        }

        def ch_existing_data = Channel.fromPath("${param_path}/${file_pattern}")
            .map { file ->
                def sample_id = file.baseName.replaceAll("\\Q${file_suffix}\\E\$", '')
                [ [id: sample_id], file ]
            }

        // Apply transformation function if provided
        if ( transform_fn ) {
            ch_existing_data = ch_existing_data.map(transform_fn)
        }

        return ch_existing_data
    }
}


def createExistingFileChannel(param_path, transform_fn = null) {

    // Create existing data channel
    if ( !param_path ) {
        return Channel.empty()
    } else {
        def file_path = file(param_path)

        if ( !file_path.isFile() ) {
            error "The provided path is not a file: ${file_path}"
        }

        def ch_existing_data = Channel.fromPath( param_path, checkIfExists: true )

        // Apply transformation if provided
        if ( transform_fn ) {
            ch_existing_data = ch_existing_data.map( transform_fn )
        } 

        return ch_existing_data
    }
}