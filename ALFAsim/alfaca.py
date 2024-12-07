from pathlib import Path
from alfasim_sdk import convert_description_to_alfacase

case_description = CaseDescription(  [...] )
alfacase_content = convert_description_to_alfacase(case_description)

# Dump the content to a file
alfacase_file = Path("c:\\user\\") / 'my_project.alfacase'
alfacase_file.write_text(data=alfacase_content, encoding='UTF-8')