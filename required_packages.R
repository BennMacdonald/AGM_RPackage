required_packages <- function(required_packages)
{
	remaining_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

	if(length(remaining_packages)) 
	{
		install.packages(remaining_packages)
	}

	for(package_name in required_packages)
	{
		library(package_name,character.only=TRUE)
	}
}

