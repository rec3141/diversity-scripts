echo "sed -E 's#\([0-9]+\)##g' $1 > $2"
sed -E 's#\([0-9]+\)##g' $1 > $2
head $2
echo "
--
done"
