package br.unisinos;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class Final {
	static void bubbleSort(int arr[], int n) {
		int i, j, temp;
		boolean swapped;
		for (i = 0; i < n - 1; i++) {
			swapped = false;
			for (j = 0; j < n - i - 1; j++) {
				if (arr[j] > arr[j + 1]) {

					// Swap arr[j] and arr[j+1]
					temp = arr[j];
					arr[j] = arr[j + 1];
					arr[j + 1] = temp;
					swapped = true;
				}
			}

			if (swapped == false)
				break;
		}
	}

	static void insertionsort(int arr[]) {
		int n = arr.length;
		for (int i = 1; i < n; ++i) {
			int key = arr[i];
			int j = i - 1;

			while (j >= 0 && arr[j] > key) {
				arr[j + 1] = arr[j];
				j = j - 1;
			}
			arr[j + 1] = key;
		}
	}

	static void selectionSort(int arr[]) {
		int n = arr.length;
		for (int i = 0; i < n - 1; i++) {
			int min_idx = i;
			for (int j = i + 1; j < n; j++)
				if (arr[j] < arr[min_idx])
					min_idx = j;
			int temp = arr[min_idx];
			arr[min_idx] = arr[i];
			arr[i] = temp;
		}
	}

	static void heapSort(int arr[]) {
		int N = arr.length;
		for (int i = N / 2 - 1; i >= 0; i--)
			heapify(arr, N, i);
		for (int i = N - 1; i > 0; i--) {
			int temp = arr[0];
			arr[0] = arr[i];
			arr[i] = temp;
			heapify(arr, i, 0);
		}
	}

	static void heapify(int arr[], int N, int i) {
		int largest = i;
		int l = 2 * i + 1;
		int r = 2 * i + 2;
		if (l < N && arr[l] > arr[largest])
			largest = l;
		if (r < N && arr[r] > arr[largest])
			largest = r;
		if (largest != i) {
			int swap = arr[i];
			arr[i] = arr[largest];
			arr[largest] = swap;
			heapify(arr, N, largest);
		}
	}

	static int shellSort(int arr[]) {
		int n = arr.length;
		for (int gap = n / 2; gap > 0; gap /= 2) {
			for (int i = gap; i < n; i += 1) {
				int temp = arr[i];
				int j;
				for (j = i; j >= gap && arr[j - gap] > temp; j -= gap)
					arr[j] = arr[j - gap];
				arr[j] = temp;
			}
		}
		return 0;
	}

	static void merge(int arr[], int l, int m, int r) {
		int n1 = m - l + 1;
		int n2 = r - m;
		int L[] = new int[n1];
		int R[] = new int[n2];
		for (int i = 0; i < n1; ++i)
			L[i] = arr[l + i];
		for (int j = 0; j < n2; ++j)
			R[j] = arr[m + 1 + j];
		int i = 0, j = 0;
		int k = l;
		while (i < n1 && j < n2) {
			if (L[i] <= R[j]) {
				arr[k] = L[i];
				i++;
			} else {
				arr[k] = R[j];
				j++;
			}
			k++;
		}
		while (i < n1) {
			arr[k] = L[i];
			i++;
			k++;
		}
		while (j < n2) {
			arr[k] = R[j];
			j++;
			k++;
		}
	}

	static void mergeSort(int arr[], int l, int r) {
		if (l < r) {
			int m = l + (r - l) / 2;
			mergeSort(arr, l, m);
			mergeSort(arr, m + 1, r);
			merge(arr, l, m, r);
		}
	}

	// QuickSort
	static void swap(int[] arr, int i, int j) {
		int temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
	}

	static int partition(int[] arr, int low, int high) {
		int pivot = arr[high];
		int i = (low - 1);

		for (int j = low; j <= high - 1; j++) {
			if (arr[j] < pivot) {
				i++;
				swap(arr, i, j);
			}
		}
		swap(arr, i + 1, high);
		return (i + 1);
	}

	static void quickSort(int[] arr, int low, int high) {
		if (low < high) {
			int pi = partition(arr, low, high);
			quickSort(arr, low, pi - 1);
			quickSort(arr, pi + 1, high);
		}
	}

	public static double calculateMean(long[] times) {
		double sum = 0.0;
		for (long time : times) {
			sum += time;
		}
		return sum / times.length;
	}

	public static double calculateVariance(long[] times, double mean) {
		double sum = 0.0;
		for (long time : times) {
			sum += Math.pow(time - mean, 2);
		}
		return sum / times.length;
	}

	public static double calculateStandardDeviation(double variance) {
		return Math.sqrt(variance);
	}

	public static double calculateFilteredMean(long[] times, double mean, double stdDev) {
		List<Long> filteredTimes = new ArrayList<>();
		for (long time : times) {
			if (Math.abs(time - mean) <= stdDev) {
				filteredTimes.add(time);
			}
		}
		long sum = 0;
		for (long time : filteredTimes) {
			sum += time;
		}
		return (double) sum / filteredTimes.size();
	}

	public static void main(String[] args) {
		int arraySize = 65536;
		int[] arrayEmOrdemCrescenteSemRepetidos = new int[arraySize];
		for (int i = 0; i < arrayEmOrdemCrescenteSemRepetidos.length; i++) {
			arrayEmOrdemCrescenteSemRepetidos[i] = i + 1;
		}

		int[] arrayEmOrdemDecrescenteSemRepetidos = new int[arraySize];
		for (int i = 0; i < arrayEmOrdemDecrescenteSemRepetidos.length; i++) {
			arrayEmOrdemDecrescenteSemRepetidos[i] = arraySize - i;
		}

		List<Integer> numerosEmbaralhados = new ArrayList<>();
		for (int i = 1; i <= arraySize; i++) {
			numerosEmbaralhados.add(i);
		}
		Collections.shuffle(numerosEmbaralhados);
		int[] arrayAleatorioSemRepetidos = new int[arraySize];
		for (int i = 0; i < arrayAleatorioSemRepetidos.length; i++) {
			arrayAleatorioSemRepetidos[i] = numerosEmbaralhados.get(i);
		}

		Random booleano = new Random();
		List<Integer> numerosEmbaralhadosRepetidos = new ArrayList<>();
		for (int i = 1; i <= arraySize; i++) {
			numerosEmbaralhadosRepetidos.add(booleano.nextInt(arraySize) + 1);
		}
		Collections.shuffle(numerosEmbaralhadosRepetidos);
		int[] arrayAleatorioComRepetidos = new int[arraySize];
		for (int i = 0; i < arrayAleatorioComRepetidos.length; i++) {
			arrayAleatorioComRepetidos[i] = numerosEmbaralhadosRepetidos.get(i);
		}

		int[] arrayParaBubbleSort = Arrays.copyOf(arrayEmOrdemCrescenteSemRepetidos, arraySize);
		int[] arrayParaInsertionSort = Arrays.copyOf(arrayEmOrdemCrescenteSemRepetidos, arraySize);
		int[] arrayParaSelectionSort = Arrays.copyOf(arrayEmOrdemCrescenteSemRepetidos, arraySize);
		int[] arrayParaHeapSort = Arrays.copyOf(arrayEmOrdemCrescenteSemRepetidos, arraySize);
		int[] arrayParaShellSort = Arrays.copyOf(arrayEmOrdemCrescenteSemRepetidos, arraySize);
		int[] arrayParaMergeSort = Arrays.copyOf(arrayEmOrdemCrescenteSemRepetidos, arraySize);
		int[] arrayParaQuickSort = Arrays.copyOf(arrayEmOrdemCrescenteSemRepetidos, arraySize);
		
		DecimalFormat df = new DecimalFormat("#");
        df.setMaximumFractionDigits(10);

		long startTime, endTime;
		long[] durations = new long[10];
		double mean, variance, stdDev, filteredMean;

		// Bubble Sort
		for (int i = 0; i < 10; i++) {
			int[] arrayCopy = Arrays.copyOf(arrayParaBubbleSort, arraySize);
			startTime = System.nanoTime();
			bubbleSort(arrayCopy, arrayCopy.length);
			endTime = System.nanoTime();
			durations[i] = endTime - startTime;
		}
		mean = calculateMean(durations);
		variance = calculateVariance(durations, mean);
		stdDev = calculateStandardDeviation(variance);
		filteredMean = calculateFilteredMean(durations, mean, stdDev);
		System.out.println("Bubble Sort - Média Filtrada (nanosegundos): " + df.format(filteredMean));

		// Insertion Sort
		for (int i = 0; i < 10; i++) {
			int[] arrayCopy = Arrays.copyOf(arrayParaInsertionSort, arraySize);
			startTime = System.nanoTime();
			insertionsort(arrayCopy);
			endTime = System.nanoTime();
			durations[i] = endTime - startTime;
		}
		mean = calculateMean(durations);
		variance = calculateVariance(durations, mean);
		stdDev = calculateStandardDeviation(variance);
		filteredMean = calculateFilteredMean(durations, mean, stdDev);
		System.out.println("Insertion Sort - Média Filtrada (nanosegundos): " + df.format(filteredMean));

		// Selection Sort
		for (int i = 0; i < 10; i++) {
			int[] arrayCopy = Arrays.copyOf(arrayParaSelectionSort, arraySize);
			startTime = System.nanoTime();
			selectionSort(arrayCopy);
			endTime = System.nanoTime();
			durations[i] = endTime - startTime;
		}
		mean = calculateMean(durations);
		variance = calculateVariance(durations, mean);
		stdDev = calculateStandardDeviation(variance);
		filteredMean = calculateFilteredMean(durations, mean, stdDev);
		System.out.println("Selection Sort - Média Filtrada (nanosegundos): " + df.format(filteredMean));

		// Heap Sort
		for (int i = 0; i < 10; i++) {
			int[] arrayCopy = Arrays.copyOf(arrayParaHeapSort, arraySize);
			startTime = System.nanoTime();
			heapSort(arrayCopy);
			endTime = System.nanoTime();
			durations[i] = endTime - startTime;
		}
		mean = calculateMean(durations);
		variance = calculateVariance(durations, mean);
		stdDev = calculateStandardDeviation(variance);
		filteredMean = calculateFilteredMean(durations, mean, stdDev);
		System.out.println("Heap Sort - Média Filtrada (nanosegundos): " + df.format(filteredMean));

		// Shell Sort
		for (int i = 0; i < 10; i++) {
			int[] arrayCopy = Arrays.copyOf(arrayParaShellSort, arraySize);
			startTime = System.nanoTime();
			shellSort(arrayCopy);
			endTime = System.nanoTime();
			durations[i] = endTime - startTime;
		}
		mean = calculateMean(durations);
		variance = calculateVariance(durations, mean);
		stdDev = calculateStandardDeviation(variance);
		filteredMean = calculateFilteredMean(durations, mean, stdDev);
		System.out.println("Shell Sort - Média Filtrada (nanosegundos): " + df.format(filteredMean));

		// Merge Sort
		for (int i = 0; i < 10; i++) {
			int[] arrayCopy = Arrays.copyOf(arrayParaMergeSort, arraySize);
			startTime = System.nanoTime();
			mergeSort(arrayCopy, 0, arrayCopy.length - 1);
			endTime = System.nanoTime();
			durations[i] = endTime - startTime;
		}
		mean = calculateMean(durations);
		variance = calculateVariance(durations, mean);
		stdDev = calculateStandardDeviation(variance);
		filteredMean = calculateFilteredMean(durations, mean, stdDev);
		System.out.println("Merge Sort - Média Filtrada (nanosegundos): " + df.format(filteredMean));

		// Quick Sort
		for (int i = 0; i < 10; i++) {
			int[] arrayCopy = Arrays.copyOf(arrayParaQuickSort, arraySize);
			startTime = System.nanoTime();
			quickSort(arrayCopy, 0, arrayCopy.length - 1);
			endTime = System.nanoTime();
			durations[i] = endTime - startTime;
		}
		mean = calculateMean(durations);
		variance = calculateVariance(durations, mean);
		stdDev = calculateStandardDeviation(variance);
		filteredMean = calculateFilteredMean(durations, mean, stdDev);
		System.out.println("Quick Sort - Média Filtrada (nanosegundos): " + df.format(filteredMean));
	}
}
