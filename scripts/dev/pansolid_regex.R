
pansolid_filename_regex <- regex(
  
  r"[
  ^Annotated            # Standard beginning
  (_|                   # Option if no panel in filename
  _v2.+_)               # Option if panel name included
  (WS\d{6})             # Worksheet
  _
  (\d{6,8})             # Lab number (usually 8 digits but rare cases can be 6 digits)
  (|a|b|c|d)            # Letter suffixes used for repeat testing
  _
  ([:alnum:]{3,30})     # Patient name - alphanumeric characters only
  (.xlsx|
  _S.+.xlsx|
  _S.+|
  _CNV_processed.xlsx)   # Ending varies between patients and controls
  ]",
  
  comments = TRUE
  
)

pansolid_filepath_regex <- regex(
  
  r"[
  .+
  WorksheetAnalysedData/
  (WS\d{6})                    # Worksheet folder
  /                             
  (\w{1,30})                   # Panel folder
  /.+
  ]",
  
  comments = TRUE
  
)




