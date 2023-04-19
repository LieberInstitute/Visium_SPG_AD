library('here')

aws_url = ''
project_name = 'Visium_SPG_AD'
samui_url = 'https://samuibrowser.com'

#   Get Samui sample IDs simply from existing files rather than reading in and
#   cleaning the sample-info excel sheet
samui_ids = list.files(
    here('processed-data', '16_samui'),
    pattern = '^V.*_[A-D]1_Br[0-9]{4}_(AD|control)$'
)

#   We want to show the Br3880 D1 sample by default when loading Samui (so we
#   put it first) 
is_first_id = grepl('Br3880_D1', samui_ids)
samui_ids = c(samui_ids[is_first_id], samui_ids[!is_first_id])

#   Construct the full link we'll share
full_url = paste0(
    samui_url, '/from?url=', aws_url, '/', project_name, '/&s=',
    paste(samui_ids, sep = '&s=')
)

print('Full URL:')
print(full_url)
