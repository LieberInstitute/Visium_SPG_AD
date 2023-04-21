library('here')

aws_url = '[AWS URL HERE]'
project_name = 'Visium_SPG_AD'
samui_url = 'https://samuibrowser.com'
first_id = 'V10A27106_D1_Br3880_AD'

#   Get Samui sample IDs simply from existing files rather than reading in and
#   cleaning the sample-info excel sheet
samui_ids = list.files(
    here('processed-data', '16_samui'),
    pattern = '^V.*_[A-D]1_Br[0-9]{4}_(AD|control)$'
)

stopifnot(first_id %in% samui_ids)
samui_ids = c(first_id, samui_ids[samui_ids != first_id])

#   Construct the full link we'll share
full_url = paste0(
    samui_url, '/from?url=', aws_url, '/', project_name, '/&s=',
    paste(samui_ids, collapse = '&s=')
)

print('Full URL:')
print(full_url)
