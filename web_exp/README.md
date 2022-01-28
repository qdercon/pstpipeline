# Web experiment

### includes demographics, allocation, PST, and questionnaires

The web experiment can be run in [JATOS](https://www.jatos.org/) 
simply by importing the compressed folder. Running outside of
JATOS will require some modification to the scripts, as will running locally outside of Prolific due to automatic recording and redirects for Prolific IDs (these changes are minor and are commented out in the scripts).

In addition, using the allocation script will require some extra info to be added in JATOS. It is adapted from this example study: https://github.com/JATOS/JATOS_examples/raw/master/examples/randomize_tasks_between_workers.jzip

First, create a new study batch under 'Worker & Batch Manager' in the JATOS GUI (or just use the 'Default' batch). Then, click on 'Properties', and add an array like the following to the 'JSON input': 

```
{
  "conditionCounts": {
    "A": 20,
    "B": 20
  },
  "completion_link": "https://www.mrc-cbu.cam.ac.uk/"
}
```

The first time the allocation component is run, an array containing 20 'A's and 20 'B's will be added to the Batch Session Data, and one will be randomly selected - A will result in the non-distanced task being run, and B the distanced. Importantly, whichever letter was selected from the array will then be removed. Subsequently, other workers will have their allocation taken from the array (the script checks if one already exists). Once all 40 workers have been run, the Batch Session Data will need to be cleared, otherwise no more workers will be able to run. 

The completion link is required for Prolific to mark the experiment as complete, and is used by the final questionnaire questions script to redirect after all results are logged. This can be an arbitrary link as seen above.

Please note that all credit for the digit span task plugin/script goes to: https://github.com/mahiluthra/working_memory_tests

