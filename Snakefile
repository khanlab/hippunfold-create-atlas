from snakebids import bids

configfile: 'config.yml'


subjects = config['subjects']

rule all:
    input:
         label = expand(bids(
                den="{density}",
                suffix="maxprob.label.gii",
                desc='freesurfersubfields',
                space="{space}",
                label="hipp",
                hemi='{hemi}',
                subject='all'),
                density='unfoldiso',
                space='corobl',
                hemi=['L','R'])

hemi_lookup={'L': 'lh','R': 'rh'}


rule import_subfields:
    input:
        lambda wildcards: config['fs_subfields'].format(subject=wildcards.subject, hemi=hemi_lookup[wildcards.hemi])
    output:
        bids(
            desc='freesurfersubfields',
            suffix="dseg.nii.gz",
            hemi='{hemi}',
            subject='{subject}'
        ),
    shell: 
        'mri_convert {input} {output}'


rule convert_fs_xfm:
    input:
        config['fs_xfm_t1_t2']
    output:        
        t1_to_t2 = bids(
            from_='T1w',
            to_='T2w',
            suffix="xfm.txt",
            type_='itk',
            subject='{subject}'
        ),
        t2_to_t1 = bids(
            from_='T2w',
            to_='T1w',
            suffix="xfm.txt",
            type_='itk',
            subject='{subject}'
        ),

    shell:
        'lta_convert --inlta {input} --outitk {output.t1_to_t2} && '
        'lta_convert --inlta {input}  --outitk {output.t2_to_t1} --invert '


#not needed if we just use surfaces in T2w space??
rule warp_seg_to_corobl_crop:
    input:
        nii=rules.import_subfields.output,
        xfm_t2_corobl=bids(
            root=config['in_hippunfold']+'/work',
            datatype="warps",
            subject='{subject}',
            suffix="xfm.txt",
            from_="T2w",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
        ref=bids(
                root=config['in_hippunfold']+'/work',
                datatype='anat',
                subject='{subject}',
                hemi='{hemi}',
                space='corobl',
                desc='preproc',
                suffix='T2w.nii.gz'),
        xfm_t1_t2 = bids(
            from_='T1w',
            to_='T2w',
            suffix="xfm.txt",
            type_='itk',
            subject='{subject}'
        ),


    output:
        nii=bids(
            desc='freesurfersubfields',
            suffix="dseg.nii.gz",
            space='corobl',
            hemi='{hemi}',
            subject='{subject}')
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.ref}  -t {input.xfm_t2_corobl} {input.xfm_t1_t2}"


#convert to workbench label volume (nifti with label info stored in header)
rule dseg_to_label_vol:
    input:
        nii=rules.warp_seg_to_corobl_crop.output.nii,
        label_list=config['subfield_label_list']
    output:
        nii=bids(
            desc='freesurfersubfields',
            suffix="dsegwb.nii.gz",
            space='{space}',
            hemi='{hemi}',
            subject='{subject}')
    shell:
         "wb_command -volume-label-import {input.nii} {input.label_list} {output.nii}"


rule map_subfields_to_surf:
    input:
        nii=rules.dseg_to_label_vol.output.nii,
        surf=bids(
            root=config['in_hippunfold']+'/work',
            datatype="surf",
            den="{density}",
            suffix="midthickness.surf.gii",
            space="{space}",
            label="hipp",
            hemi="{hemi}",
            subject='{subject}')
    output:
        label_gii=bids(
            den="{density}",
            suffix="subfields.label.gii",
            desc='freesurfersubfields',
            space="{space}",
            label="hipp",
            hemi="{hemi}",
            subject='{subject}')
    shell:
        "wb_command -volume-label-to-surface-mapping {input.nii} {input.surf} {output.label_gii}"

rule create_cifti:
    input:
        left_label=bids(
            den="{density}",
            suffix="subfields.label.gii",
            desc='freesurfersubfields',
            space="{space}",
            label="hipp",
            hemi="L",
            subject='{subject}'
        ),
        right_label=bids(
            den="{density}",
            suffix="subfields.label.gii",
            desc='freesurfersubfields',
            space="{space}",
            label="hipp",
            hemi="R",
            subject='{subject}'
        ),
    output:
        cifti=bids(
            den="{density}",
            suffix="dseg.dlabel.nii",
            desc='freesurfersubfields',
            space="{space}",
            label="hipp",
            subject='{subject}'
        ),
    group:
        "subj"
    shell:
        "wb_command  -cifti-create-label {output}"
        " -left-label {input.left_label}"
        " -right-label {input.right_label}"

rule merge_cifti:
    input:
        ciftis = expand(bids(
                den="{density}",
                suffix="dseg.dlabel.nii",
                desc='freesurfersubfields',
                space="{space}",
                label="hipp",
                subject='{subject}'),
            subject=subjects,
            allow_missing=True,
            ),
    output:
        cifti = bids(
                den="{density}",
                suffix="merged.dlabel.nii",
                desc='freesurfersubfields',
                space="{space}",
                label="hipp",
                subject='all')
    params:
        in_ciftis = lambda wildcards, input: ' '.join( [ f'-cifti {i}' for i in input])
    shell:
        'wb_command -cifti-merge {output} {params.in_ciftis}'

rule create_prob_label:
    input:
        cifti = rules.merge_cifti.output
    output:
        cifti = bids(
                den="{density}",
                suffix="probseg.dscalar.nii",
                desc='freesurfersubfields',
                space="{space}",
                label="hipp",
                subject='all')
    shell:
        'wb_command -cifti-label-probability {input} {output} -exclude-unlabeled'

rule create_maxprob:
    input: 
        cifti = rules.create_prob_label.output.cifti,
        label_list=config['subfield_label_list']
    output:
        cifti_dscalar = temp(bids(
                den="{density}",
                suffix="maxprob.dscalar.nii",
                desc='freesurfersubfields',
                space="{space}",
                label="hipp",
                subject='all')),
        cifti_dlabel = bids(
                den="{density}",
                suffix="maxprob.dlabel.nii",
                desc='freesurfersubfields',
                space="{space}",
                label="hipp",
                subject='all')

    shell:
        'wb_command -cifti-reduce {input.cifti} INDEXMAX {output.cifti_dscalar} && '
        'wb_command -cifti-label-import {output.cifti_dscalar} {input.label_list} {output.cifti_dlabel}'


rule cifti_to_gifti:
    """ to allow resampling to diff spaces """
    input:
        cifti_dlabel = bids(
                den="{density}",
                suffix="maxprob.dlabel.nii",
                desc='freesurfersubfields',
                space="{space}",
                label="hipp",
                subject='all')
    output:
        label_lh = bids(
                den="{density}",
                suffix="maxprob.label.gii",
                desc='freesurfersubfields',
                space="{space}",
                label="hipp",
                hemi='L',
                subject='all'),
        label_rh = bids(
                den="{density}",
                suffix="maxprob.label.gii",
                desc='freesurfersubfields',
                space="{space}",
                label="hipp",
                hemi='R',
                subject='all'),
    shell:
        'wb_command -cifti-separate {input} COLUMN '
        ' -label CORTEX_LEFT {output.label_lh} '
        ' -label CORTEX_RIGHT {output.label_rh} '


rule resample_to_unfoldiso:
    input:
        label = bids(
                den="0p5mm",
                suffix="maxprob.label.gii",
                desc='freesurfersubfields',
                space="{space}",
                label="hipp",
                hemi='{hemi}',
                subject='all'),
        current = config['ref_surf_unfold'].format(density='0p5mm'),
        new = lambda wildcards: config['ref_surf_unfold'].format(density=wildcards.density),
    params:
        method = 'BARYCENTRIC'
    output:
        label = bids(
                den="{density,unfoldiso}",
                suffix="maxprob.label.gii",
                desc='freesurfersubfields',
                space="{space}",
                label="hipp",
                hemi='{hemi}',
                subject='all'),
    shell:
        "wb_command -label-resample {input.label} {input.current} {input.new} {params.method} {output.label} -bypass-sphere-check"


        
