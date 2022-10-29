#include<stdio.h>
 
// Function implementing the Counting Sort
void counting_sort(int arr[], int n){
 
    // Step 1: Finding the maximum element
    int maximum = -1;
 
    for(int i=0;i<n;i++){
        if(maximum < arr[i])
            maximum = arr[i];
    }
 
    // Step 2: Initialize a count array of length = (max + 1) 
 
    int count[maximum+1];
    for(int i=0;i<=maximum;i++){
        count[i] = 0;
    }
 
    // Step 3: Fill the count array accordingly
 
    for(int i=0;i<n;i++){
        count[arr[i]]++;
    }
 
    // Step 4: Calculate cumulative frequency in Count Array
 
    for(int i=1;i<=maximum;i++){
        printf("%d\n", count[i-1]);
        count[i] += count[i-1];
    }
 
    // Step 5: Fix the values in the sorted array
 
    int sorted_arr[n];
 
    for(int i=0;i<n;i++){
        sorted_arr[count[arr[i]]-1] = arr[i];
        count[arr[i]]--;
    }
 
    // Printing the sorted array
 
    for(int i=0;i<n;i++){
        printf("%d (c: %d)\n", sorted_arr[i], count[i]);
    }
    printf("\n");
}
 
// The main function
int main(){
 
    // Initializing the array
    int arr[] = { 5, 2, 5, 3, 6, 1, 5, 3, 9, 6};
 
    // Size of the array
    int n = sizeof(arr)/sizeof(arr[0]);
 
    // Function call to the sorting function
    counting_sort(arr, n);
    return 1;
}