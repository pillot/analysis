#!/bin/bash

ARGS_ALL="--session ${OVERRIDE_SESSION:-default} --shm-segment-size $SHMSIZE"
if [[ $NTIMEFRAMES == -1 ]]; then NTIMEFRAMES_CMD= ; else NTIMEFRAMES_CMD="--max-tf $NTIMEFRAMES"; fi
DISABLE_DIGIT_ROOT_INPUT="--disable-root-input"
if [[ -z "$CONFIG_EXTRA_PROCESS_o2_mch_reco_workflow" ]]; then MCH_CONFIG= ; else MCH_CONFIG="--configKeyValues \"$CONFIG_EXTRA_PROCESS_o2_mch_reco_workflow\""; fi

WORKFLOW="o2-ctf-reader-workflow $ARGS_ALL --delay $TFDELAY $NTIMEFRAMES_CMD --ctf-input $INPUT_FILE_LIST ${INPUT_FILE_COPY_CMD+--copy-cmd} ${INPUT_FILE_COPY_CMD} --ctf-dict ctf_dictionary.root --onlyDet $WORKFLOW_DETECTORS $ARGS_EXTRA_PROCESS_o2_ctf_reader_workflow"
[[ -z "$DISABLE_ROOT_OUTPUT" ]] && WORKFLOW+=" | o2-tfidinfo-writer-workflow $ARGS_ALL"
WORKFLOW+=" | o2-mid-reco-workflow $ARGS_ALL $DISABLE_ROOT_OUTPUT $DISABLE_MC"
WORKFLOW+=" | o2-mch-reco-workflow $ARGS_ALL $DISABLE_DIGIT_ROOT_INPUT $DISABLE_ROOT_OUTPUT $DISABLE_MC $ARGS_EXTRA_PROCESS_o2_mch_reco_workflow $MCH_CONFIG"
WORKFLOW+=" | o2-muon-tracks-matcher-workflow $ARGS_ALL"
WORKFLOW+=" | o2-muon-tracks-writer-workflow $ARGS_ALL"
WORKFLOW+=" | o2-dpl-run $ARGS_ALL --run"

[[ $WORKFLOWMODE == "print" ]] && echo "#Workflow command:\n\n${WORKFLOW}\n" | sed -e "s/\\\\n/\n/g" -e"s/| */| \\\\\n/g" | eval cat $( [[ $WORKFLOWMODE == "dds" ]] && echo '1>&2')
[[ $WORKFLOWMODE != "print" ]] && eval $WORKFLOW
