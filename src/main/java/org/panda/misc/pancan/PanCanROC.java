package org.panda.misc.pancan;

import org.panda.misc.MutexReader;
import org.panda.misc.MutexReader.*;
import org.panda.resource.CancerGeneBushman;
import org.panda.resource.CancerGeneCensus;
import org.panda.resource.GO;
import org.panda.resource.OncoKB;
import org.panda.utility.CollectionUtil;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class PanCanROC
{
	static String testDir = "/home/babur/Documents/PanCan/PanCan-results";
	static String ctrlDir = "/home/babur/Documents/PanCan/PanCan-shuffled-?-results";

	public static void main(String[] args) throws IOException
	{
		doReachAssessment();
//		findBestMaxDiv();
	}

	private static void doReachAssessment()
	{
		int div = 30;

		Map<String, MutexReader.DirectoryFilter> filterMap = new LinkedHashMap<>();
		filterMap.put("W", f -> !hasNameOverThan(f.getName(), div) && !f.getName().endsWith("PC2v8"));
		filterMap.put("P", f -> !hasNameOverThan(f.getName(), div) && !(f.getName().equals("no-network") || f.getName().startsWith("REACH") || f.getName().startsWith("r")));
		filterMap.put("R", f -> !hasNameOverThan(f.getName(), div) && !(f.getName().equals("no-network") || f.getName().equals("PC2v8") || f.getName().startsWith("r")));
		filterMap.put("WP", f -> !hasNameOverThan(f.getName(), div) && !(f.getName().startsWith("REACH") || f.getName().startsWith("r")));
		filterMap.put("WR", f -> !hasNameOverThan(f.getName(), div) && !(f.getName().startsWith("PC") || f.getName().startsWith("r")));
		filterMap.put("PR", f -> !hasNameOverThan(f.getName(), div) && !(f.getName().equals("no-network") || f.getName().startsWith("r")));
		filterMap.put("WPR", f -> !hasNameOverThan(f.getName(), div) && ! f.getName().startsWith("r"));
		filterMap.put("WR'", f -> !hasNameOverThan(f.getName(), div) && !(f.getName().startsWith("REACH") || f.getName().startsWith("PC")));
		filterMap.put("R'", f -> !hasNameOverThan(f.getName(), div) && !(f.getName().startsWith("REACH") || f.getName().equals("PC2v8") || f.getName().equals("no-network")));

		drawRoc(filterMap, 1000, true);
	}

	private static void findBestMaxDiv()
	{
		Map<String, MutexReader.DirectoryFilter> filterMap = new LinkedHashMap<>();
		for (int i = 1; i <= 30; i++)
		{
			int div = i;
//			filterMap.put("" + div, f -> !hasNameOverThan(f.getName(), div) && !(f.getName().startsWith("random")));
			filterMap.put("" + div, f -> (f.getPath().contains("/" + div + "/") || f.getPath().endsWith("/" + div)) && !(f.getName().startsWith("random")));
		}

		drawRoc(filterMap, 1000, true);
	}

	private static void drawRoc(Map<String, MutexReader.DirectoryFilter> filterMap, int resultSizeLimit, boolean mutex)
	{
		Map<String, List<Double>> testScoreMap = new LinkedHashMap<>();
		Map<String, List<Double>> ctrlScoreMap = new LinkedHashMap<>();

		double ss = 0;
		for (String runType : filterMap.keySet())
		{
			Object[] o = PanCanResultLoader.readGroupsWithFlattenedControl(mutex, testDir, ctrlDir, filterMap.get(runType));
			Set<Group> testResults = (Set<Group>) o[0];
			Map<String, Double> testScores = MutexReader.convertGroupsToGeneBestScores(testResults);
			List<Double> ctrlScores = (List<Double>) o[1];
			ss = (int) o[2];
			List<Double> testScoresList = new ArrayList<>(testScores.values());
			Collections.sort(testScoresList);
			Collections.sort(ctrlScores);
			testScoreMap.put(runType, testScoresList);
			ctrlScoreMap.put(runType, ctrlScores);
		}

		Map<String, Integer> cIndMap = new HashMap<>();
		Map<String, Double> minTrueMap = new HashMap<>();
		for (String type : filterMap.keySet())
		{
			cIndMap.put(type, -1);
			minTrueMap.put(type, 0D);
		}

		for (String type : filterMap.keySet())
		{
			System.out.print("\t" + type + "\t");
		}

		for (int j = 1; j <= resultSizeLimit; j++)
		{
			System.out.println();

			for (String type : filterMap.keySet())
			{
				if (testScoreMap.get(type).size() < j+1)
				{
					System.out.print("\t\t");
					continue;
				}

				while (cIndMap.get(type) < ctrlScoreMap.get(type).size() - 1 &&
					ctrlScoreMap.get(type).get(cIndMap.get(type) + 1) < testScoreMap.get(type).get(j-1))
				{
					cIndMap.put(type, cIndMap.get(type) + 1);
				}

				if (testScoreMap.get(type).get(j-1).equals(testScoreMap.get(type).get(j)))
				{
					System.out.print("\t\t");
				}
				else
				{
					double fp = (cIndMap.get(type) + 1) / ss;

					if (fp <= j)
					{
						double tp = j - fp;

						if (minTrueMap.get(type) > tp)
						{
							tp = minTrueMap.get(type);
							fp = j - tp;
						}
						else minTrueMap.put(type, tp);

						System.out.print(fp + "\t" + tp + "\t");
					}
					else System.out.print("\t\t");
				}
			}
		}
	}
	public static boolean hasNameOverThan(String name, int thr)
	{
		int x;
		try
		{
			x = Integer.parseInt(name);
		}
		catch (NumberFormatException e){return false;}

		return x > thr;
	}
}
