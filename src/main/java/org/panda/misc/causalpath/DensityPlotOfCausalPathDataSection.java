package org.panda.misc.causalpath;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Progress;
import org.panda.utility.statistics.Correlation;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.KernelDensityPlot;
import org.panda.utility.statistics.Summary;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class DensityPlotOfCausalPathDataSection
{
	public static void main(String[] args) throws IOException
	{
//		generateForCompBased("/Users/ozgun/Documents/Analyses/CausalPath-data/CPTAC-BRCA/subtypes/LumAB-vs-Basal");
//		printMinimumChangeForCompBased();
//		printMinimumChangeForCompBasedPaired();
		generateForCorrBased("/Users/ozgun/Documents/Analyses/CausalPath-data/CPTAC-BRCA/correlation-based-phospho");
//		findProblematicDataFiles();
//		findSmallestCorrelation();
	}

	private static void generateForCompBased(String dir) throws IOException
	{
		Set<String> selected = getSelectedIDsFromResults(dir + "/results.txt");

		String dataFile = dir + "/" + FileUtil.lines(dir + "/parameters.txt")
			.filter(l -> l.startsWith("proteomics-values-file ="))
			.map(l -> l.split("=")[1].trim()).findFirst().get();

		Set<String> ctrl = FileUtil.lines(dir + "/parameters.txt")
			.filter(l -> l.startsWith("control-value-column ="))
			.map(l -> l.split("=")[1].trim()).collect(Collectors.toSet());
		Set<String> test = FileUtil.lines(dir + "/parameters.txt")
			.filter(l -> l.startsWith("test-value-column ="))
			.map(l -> l.split("=")[1].trim()).collect(Collectors.toSet());

		Map<String, double[][]> data = readDataForComp(dataFile, ctrl, test);
		Map<String, double[]> allData = readAllData(dataFile);

		double min = Double.MAX_VALUE;
		Map<String, Double> valMap = new HashMap<>();
		for (String id : data.keySet())
		{
			double[] vc = data.get(id)[0];
			double[] vt = data.get(id)[1];

			double[] all = allData.get(id);

			double diff = Summary.mean(vt) - Summary.mean(vc);
			diff /= Summary.stdev(all);
			if (!Double.isNaN(diff))
			{
				valMap.put(id, diff);
				if (selected.contains(id) && Math.abs(diff) < min) min = Math.abs(diff);
			}
		}

		Set<String> subsetIDs = CollectionUtil.getIntersection(selected, valMap.keySet());
		printSelectionHistogram(valMap, subsetIDs);
		System.out.println("min = " + min);
	}

	private static void printMinimumChangeForCompBased() throws IOException
	{
		FileUtil.processDirsRecursive(new File("/Users/ozgun/Documents/Analyses/CausalPath-data/"), dir ->
		{
			String resultFile = dir + "/results.txt";
			if (Files.exists(Paths.get(resultFile)))
			{
				if (FileUtil.lines(dir + "/parameters.txt")	.anyMatch(l -> l.endsWith("-paired"))) return;

				String[] header = FileUtil.lines(resultFile).findFirst().get().split("\t");
				int cInd = ArrayUtil.indexOf(header, "Correlation");
				if (cInd < 0)
				{
					Set<String> selected = getSelectedIDsFromResults(resultFile);

					String dataFile = dir + "/" + FileUtil.lines(dir + "/parameters.txt")
						.filter(l -> l.startsWith("proteomics-values-file ="))
						.map(l -> l.split("=")[1].trim()).findFirst().get();

					Set<String> ctrl = FileUtil.lines(dir + "/parameters.txt")
						.filter(l -> l.startsWith("control-value-column ="))
						.map(l -> l.split("=")[1].trim()).collect(Collectors.toSet());
					Set<String> test = FileUtil.lines(dir + "/parameters.txt")
						.filter(l -> l.startsWith("test-value-column ="))
						.map(l -> l.split("=")[1].trim()).collect(Collectors.toSet());

					if (ctrl.isEmpty() || test.isEmpty()) return;

					Map<String, double[][]> data = readDataForComp(dataFile, ctrl, test);
					Map<String, double[]> allData = readAllData(dataFile);

					double min = Double.MAX_VALUE;
					for (String id : selected)
					{
						if (data.containsKey(id))
						{
							double[] vc = data.get(id)[0];
							double[] vt = data.get(id)[1];

							double[] all = allData.get(id);

							double diff = Summary.mean(vt) - Summary.mean(vc);
							double stdev = Summary.stdev(all);
							diff /= stdev;
							if (!Double.isNaN(diff))
							{
								diff = Math.abs(diff);
								if (diff < min) min = diff;
							}
						}
					}
					System.out.println(dir.getPath() + "\t" + min);
				}
			}
		});
	}

	private static void printMinimumChangeForCompBasedPaired() throws IOException
	{
		FileUtil.processDirsRecursive(new File("/Users/ozgun/Documents/Analyses/CausalPath-data/"), dir ->
		{
			String resultFile = dir + "/results.txt";
			if (Files.exists(Paths.get(resultFile)))
			{
				if (!FileUtil.lines(dir + "/parameters.txt").anyMatch(l -> l.endsWith("-paired"))) return;

				String[] header = FileUtil.lines(resultFile).findFirst().get().split("\t");
				int cInd = ArrayUtil.indexOf(header, "Correlation");
				if (cInd < 0)
				{
					Set<String> selected = getSelectedIDsFromResults(resultFile);

					String dataFile = dir + "/" + FileUtil.lines(dir + "/parameters.txt")
						.filter(l -> l.startsWith("proteomics-values-file ="))
						.map(l -> l.split("=")[1].trim()).findFirst().get();

					List<String> ctrl = FileUtil.lines(dir + "/parameters.txt")
						.filter(l -> l.startsWith("control-value-column ="))
						.map(l -> l.split("=")[1].trim()).collect(Collectors.toList());
					List<String> test = FileUtil.lines(dir + "/parameters.txt")
						.filter(l -> l.startsWith("test-value-column ="))
						.map(l -> l.split("=")[1].trim()).collect(Collectors.toList());

					if (ctrl.isEmpty() || test.isEmpty()) return;

					Map<String, double[][]> data = readDataForComp(dataFile, ctrl, test);

					double min = Double.MAX_VALUE;
					for (String id : selected)
					{
						if (data.containsKey(id))
						{
							double[] vc = data.get(id)[0];
							double[] vt = data.get(id)[1];

							double[] d = new double[vc.length];
							for (int i = 0; i < vc.length; i++)
							{
								d[i] = vt[i] - vc[i];
							}

							double diff = Summary.mean(d);
							double stdev = Summary.stdev(d);
							diff /= stdev;
							if (!Double.isNaN(diff))
							{
								diff = Math.abs(diff);
								if (diff < min) min = diff;
							}
						}
					}
					System.out.println(dir.getPath() + "\t" + min);
				}
			}
		});
	}

	private static void printSelectionHistogram(Map<String, Double> valMap, Set<String> selection)
	{
		double[] sum = new double[]{0};
		Histogram h1 = new Histogram(0.2);
		valMap.values().forEach(v -> h1.count(v));
		valMap.values().forEach(v -> sum[0] += v);

		Histogram h2 = new Histogram(0.2);
		selection.forEach(id -> h2.count(valMap.get(id)));

		h1.printTogether(h2);

		System.out.println("sum = " + sum[0]);
		System.out.println("average = " + (sum[0] / h1.getTotal()));
	}

	private static Map<String, Map<String, Double>> loadValueChanges(String file) throws IOException
	{
		Map<String, Map<String, Double>> maps = new HashMap<>();
		Map<String, Double> current = null;
		int pInd = -1;

		Scanner sc = new Scanner(Paths.get(file));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (line.startsWith("Data type:"))
			{
				current = null;
				String type = line.substring(line.lastIndexOf(".") + 1);
				String[] header = sc.nextLine().split("\t");
				pInd = ArrayUtil.indexOf(header, "P-value");
				if (pInd > 0)
				{
					current = new HashMap<>();
					maps.put(type, current);
				}
			}
			else if (!line.isEmpty() && current !=  null && pInd > 0)
			{
				String[] t = line.split("\t");
				if (!t[pInd].equals("NaN"))
				{
					current.put(t[0], Double.valueOf(t[pInd]));
				}
			}
		}
		return maps;
	}

	private static Set<String> getSelectedIDsFromResults(String file)
	{
		Set<String> ids = new HashSet<>();

		String[] h = FileUtil.lines(file).findFirst().get().split("\t");
		int sInd = ArrayUtil.indexOf(h, "Source data ID");
		int tInd = ArrayUtil.indexOf(h, "Target data ID");
		FileUtil.lines(file).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			ids.add(t[sInd]);
			ids.add(t[tInd]);
		});
		return ids;
	}

	private static Set<String> getSelectedIDsFromResultsCorr(String file)
	{
		Set<String> ids = new HashSet<>();

		String[] h = FileUtil.lines(file).findFirst().get().split("\t");
		int sInd = ArrayUtil.indexOf(h, "Source data ID");
		int tInd = ArrayUtil.indexOf(h, "Target data ID");
		FileUtil.lines(file).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[sInd].compareTo(t[tInd]) < 0 ? t[sInd] + ":" + t[tInd] : t[tInd] + ":" + t[sInd];
			ids.add(id);
		});
		return ids;
	}

	private static void generateForCorrBased(String dir)
	{
		String dataFile = dir + "/" + FileUtil.lines(dir + "/parameters.txt")
			.filter(l -> l.startsWith("proteomics-values-file ="))
			.map(l -> l.split("=")[1].trim()).findFirst().get();

		Set<String> columns = FileUtil.lines(dir + "/parameters.txt")
			.filter(l -> l.startsWith("value-column ="))
			.map(l -> l.split("=")[1].trim()).collect(Collectors.toSet());

		List<String> columnsList = FileUtil.lines(dir + "/parameters.txt")
			.filter(l -> l.startsWith("value-column ="))
			.map(l -> l.split("=")[1].trim()).collect(Collectors.toList());

//		System.out.println("columns = " + columns.size());
//		System.out.println("columnsList = " + columnsList.size());
//		Set<String> tmp = new HashSet<>();
//		for (String s : columnsList)
//		{
//			if (tmp.contains(s)) System.out.println("s = " + s);
//			else tmp.add(s);
//		}

		Map<String, double[]> data = readDataForCorr(dataFile, columns);

		Set<String> selected = getSelectedIDsFromResultsCorr(dir + "/results.txt");

		double sumBGCorr = 0;

		double range = 0.1;
		Histogram h1 = new Histogram(range);
		Histogram h2 = new Histogram(range);

		Progress p = new Progress(data.size(), "Calculating correlations");

		for (String id1 : data.keySet())
		{
			p.tick();
			for (String id2 : data.keySet())
			{
				if (id1.compareTo(id2) < 0)
				{
					if (Math.random() > 0.01) continue;
					if (id1.split("-")[0].equals(id2.split("-")[0]))
					{
						continue;
					}
					else
					{
						System.out.print("");
					}

					String id = id1 + ":" + id2;

					double[][] v = ArrayUtil.trimNaNs(data.get(id1), data.get(id2));

					if (v[0].length >= 5)
					{
						double c = Correlation.pearsonVal(v[0], v[1]);
						if (!Double.isNaN(c))
						{
							h1.count(c);
							sumBGCorr += c;
							if (selected.contains(id)) h2.count(c);
						}
					}
				}
			}
		}

		h1.printWithSubset(h2);

		System.out.println("\nsumBGCorr = " + sumBGCorr);
		System.out.println("h1.getTotal() = " + h1.getTotal());
		System.out.println("Average correlation = " + (sumBGCorr / h1.getTotal()));
	}

	private static Map<String, double[]> readDataForCorr(String dataFile, Set<String> columns)
	{
		String[] header = FileUtil.lines(dataFile).findFirst().get().split("\t");

		Map<String, double[]> map = new HashMap<>();

		FileUtil.lines(dataFile).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[0];
			double[] v = new double[columns.size()];
			int i = 0;
			for (int j = 1; j < t.length; j++)
			{
				if (columns.contains(header[j]))
				{
					v[i++] = t[j].isEmpty() || t[j].equals("NA") ? Double.NaN : Double.valueOf(t[j]);
				}
			}
			map.put(id, v);
		});
		return map;
	}

	private static Map<String, double[]> readAllData(String dataFile)
	{
		Map<String, double[]> map = new HashMap<>();

		FileUtil.lines(dataFile).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[0];
			double[] v = new double[t.length - 4];
			for (int j = 4; j < t.length; j++)
			{
				v[j-4] = t[j].isEmpty() || t[j].equals("NA") ? Double.NaN : Double.valueOf(t[j]);
			}
			map.put(id, v);
		});
		return map;
	}

	private static Map<String, double[][]> readDataForComp(String dataFile, Set<String> ctrl, Set<String> test)
	{
		String[] header = FileUtil.lines(dataFile).findFirst().get().split("\t");

		Map<String, double[][]> map = new HashMap<>();

		FileUtil.lines(dataFile).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[0];
			double[] vc = new double[ctrl.size()];
			double[] vt = new double[test.size()];
			int i = 0;
			int j = 0;
			for (int k = 1; k < t.length; k++)
			{
				if (ctrl.contains(header[k]))
				{
					vc[i++] = t[k].isEmpty() || t[k].equals("NA") ? Double.NaN : Double.valueOf(t[k]);
				}
				if (test.contains(header[k]))
				{
					vt[j++] = t[k].isEmpty() || t[k].equals("NA") ? Double.NaN : Double.valueOf(t[k]);
				}
			}
			map.put(id, new double[][]{vc, vt});
		});
		return map;
	}

	private static Map<String, double[][]> readDataForComp(String dataFile, List<String> ctrl, List<String> test)
	{
		String[] header = FileUtil.lines(dataFile).findFirst().get().split("\t");

		Map<String, double[][]> map = new HashMap<>();

		FileUtil.lines(dataFile).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[0];
			double[] vc = new double[ctrl.size()];
			double[] vt = new double[test.size()];

			int i = 0;
			for (String colName : ctrl)
			{
				String s = t[ArrayUtil.indexOf(header, colName)];
				vc[i++] = s.isEmpty() || s.equals("NA") ? Double.NaN : Double.valueOf(s);
			}
			i = 0;
			for (String colName : test)
			{
				String s = t[ArrayUtil.indexOf(header, colName)];
				vt[i++] = s.isEmpty() || s.equals("NA") ? Double.NaN : Double.valueOf(s);
			}

			map.put(id, new double[][]{vc, vt});
		});
		return map;
	}

	private static void findProblematicDataFiles() throws IOException
	{
		FileUtil.processDirsRecursive(new File("/Users/ozgun/Documents/Analyses/CausalPath-data/"), dir ->
		{
			String paramFilename = dir.getPath() + "/parameters.txt";
			if (Files.exists(Paths.get(paramFilename)))
			{
				String dataFilename = dir + "/" + FileUtil.lines(dir + "/parameters.txt")
					.filter(l -> l.startsWith("proteomics-values-file ="))
					.map(l -> l.split("=")[1].trim()).findFirst().get();

				String[] header = FileUtil.lines(dataFilename).findFirst().get().split("\t");
				Set<String> mem = new HashSet<>();
				Set<String> rep = new HashSet<>();
				for (int i = 0; i < header.length; i++)
				{
					if (mem.contains(header[i])) rep.add(header[i]);
					else mem.add(header[i]);
				}
				if (!rep.isEmpty())
				{
					System.out.println(dir + ": " + rep);
				}
			}
		});
	}

	private static void findSmallestCorrelation() throws IOException
	{
		FileUtil.processDirsRecursive(new File("/Users/ozgun/Documents/Analyses/CausalPath-data/"), dir ->
		{
			String resultsFile = dir.getPath() + "/results.txt";
			if (Files.exists(Paths.get(resultsFile)))
			{
				String[] header = FileUtil.lines(resultsFile).findFirst().get().split("\t");
				int cInd = ArrayUtil.indexOf(header, "Correlation");
				if (cInd > 0)
				{
					Double min = FileUtil.lines(resultsFile).skip(1).map(l -> Math.abs(Double.valueOf(l.split("\t")[cInd]))).min(Double::compareTo).get();
					System.out.println(dir.getPath() + "\t" + min);
				}
			}
		});

	}
}
