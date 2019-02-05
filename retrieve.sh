while [[ $# > 1 ]]
do
	key="$1"

	case $key in
		-r|--region)
		REGION="$2"
		shift # past argument
		;;
		-d|--dates)
		DATES="$2"
		shift # past argument
		;;
		-o|--out)
		outfile="$2"
		shift # past argument
		;;
		*)
				# unknown option
		;;
	esac
	shift # past argument or value
done
python retrieve_era5.py $REGION $DATES $outfile | tee ./log/"${outfile}.log" 

