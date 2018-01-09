/* Computes the set {1,...,index_set_size}/sub_set and returns its elements in ascending order
**
**
*/
int[] index_set_diff(int[] sub_set, int sub_set_length, int index_set_size) {
  int sub_set_sorted[sub_set_length] = sort_asc(sub_set);
  int j = 1;
  int diff_set[index_set_size-sub_set_length];

  for (i in 1:index_set_size) {
    if(j<=sub_set_length && sub_set_sorted[j]==i){
      j = j + 1;
    } else {
      diff_set[i-j+1] = i;
    }
  }
  return diff_set;
}
