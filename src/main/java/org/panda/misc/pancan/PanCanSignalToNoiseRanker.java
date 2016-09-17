package org.panda.misc.pancan;

import org.panda.misc.MutexReader;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.SignalToNoiseDistance;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Ozgun Babur
 */
public class PanCanSignalToNoiseRanker
{
	public static final String BASE = "/home/babur/Documents/PanCan/";
	public static final String SIGNAL_DIR = BASE + "PanCan-results/";
	public static final String NOISE_DIR = BASE + "PanCan-shuffled-results/";

	public static void main(String[] args)
	{
		Map<String, Double> map = getDistanceMap();
		map.keySet().stream().sorted((s1, s2) -> map.get(s2).compareTo(map.get(s1))).forEach(s ->
			System.out.println(s + "\t" + map.get(s)));
	}

	private static Map<String, Double> getDistanceMap()
	{
		Map<String, Double> distMap = new HashMap<>();

		for (File level1 : new File(SIGNAL_DIR).listFiles())
		{
			if (!level1.isDirectory()) continue;
			for (File level2 : level1.listFiles())
			{
				if (!level2.isDirectory()) continue;
				for (File level3 : level2.listFiles())
				{
					if (!level3.isDirectory()) continue;
					double distance = getDistance(level3.getPath());
					String name = clipLastPart(level3.getPath());
					distMap.put(name, distance);
				}
			}
		}
		return distMap;
	}

	private static double getDistance(String dir)
	{
		String noiseDir = NOISE_DIR + clipLastPart(dir);

		Map<String, Double> signalScores = MutexReader.readBestScores(dir);
		Map<String, Double> noiseScores = MutexReader.readBestScores(noiseDir);

//		return SignalToNoiseDistance.calculate(signalScores, noiseScores);
		return getSignalStrength(signalScores, clipLastPart(dir));
	}

	private static double getSignalStrength(Map<String, Double> scores, String part)
	{
		int div = Integer.parseInt(part.substring(0, part.indexOf("/")));
		return scores.values().stream().filter(v -> v <= 0.01).count() * div;
	}

	private static String clipLastPart(String path)
	{
		int index = path.lastIndexOf("/");
		index = path.lastIndexOf("/", index-1);
		index = path.lastIndexOf("/", index-1);
		return path.substring(index + 1);
	}

}
