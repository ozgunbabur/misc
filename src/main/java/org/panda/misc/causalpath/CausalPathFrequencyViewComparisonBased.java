package org.panda.misc.causalpath;

import org.panda.resource.OncoKB;
import org.panda.resource.siteeffect.SiteEffectCollective;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class CausalPathFrequencyViewComparisonBased
{
	public static void main(String[] args) throws IOException
	{
		String base = "/Users/ozgun/Documents/Analyses/CPTAC-PanCan/";
		printGeneFeatureChanges("/Users/ozgun/Documents/Analyses/CPTAC-PanCan/clusters/against-others-diffexp", new HashSet<>(Arrays.asList("GRK1", "GRK2", "GRK3", "GRK4", "GRK5", "GRK6", "GRK7")));
//		generateSummaryGraph(base + "clusters/against-others-diffexp", path -> path.split("/")[path.split("/").length - 2], getRBFocus(), base + "sifs/RB-focus-clusters-diffexp");
	}

	private static void printGeneFeatureChanges(String base, Set<String> genes) throws IOException
	{
		Map<String, int[]> idToCnts = new HashMap<>();

		FileUtil.processDirsRecursive(new File(base), dir ->
		{
			if (!dir.getPath().contains("sitespec")) return;

			String resultFilename = dir.getPath() + "/results.txt";
			if (Files.exists(Paths.get(resultFilename)))
			{
				FileUtil.linesTabbedSkip1(resultFilename).forEach(t ->
				{
					if (genes.contains(t[0]))
					{
						String id = t[4];
						boolean up = !t[5].startsWith("-");
						if (!idToCnts.containsKey(id)) idToCnts.put(id, new int[]{0, 0});
						idToCnts.get(id)[up ? 0 : 1]++;
						System.out.println(Arrays.toString(t));
					}
					if (genes.contains(t[2]))
					{
						String id = t[7];
						boolean up = !t[8].startsWith("-");
						if (!idToCnts.containsKey(id)) idToCnts.put(id, new int[]{0, 0});
						idToCnts.get(id)[up ? 0 : 1]++;
						System.out.println(Arrays.toString(t));
					}
				});
			}
		});

		SiteEffectCollective sec = new SiteEffectCollective();

		for (String id : idToCnts.keySet())
		{
			Integer eff = sec.getEffectFromID(id);
			System.out.println(id + "\t" + eff + "\t" + Arrays.toString(idToCnts.get(id)));
//			if (eff != null && eff != 0) System.out.println("map.put(\"" + id + "\", " + (-eff) + ");");
		}
	}

	private static void generateSummaryGraph(String base, NameGetter ng, Map<String, Integer> focus, String outSIFNoExt) throws IOException
	{
		Map<String, String[]> rows = new HashMap<>();
		Map<String, Integer> direcMap = new HashMap<>();
		Map<String, Set<String>> occurMap = new HashMap<>();

		FileUtil.processDirsRecursive(new File(base), dir ->
		{
			if (!dir.getPath().contains("sitespec")) return;

			String resultFilename = dir.getPath() + "/results.txt";
			if (Files.exists(Paths.get(resultFilename)))
			{
				FileUtil.linesTabbedSkip1(resultFilename).forEach(t ->
				{
					if ((focus.keySet().contains(t[4]) && focus.get(t[4]) * Double.valueOf(t[5]) > 0) || (focus.keySet().contains(t[7]) && focus.get(t[7]) * Double.valueOf(t[8]) > 0))
					{
						String key = t[4] + " " + t[7];
						if (!rows.containsKey(key)) rows.put(key, t);
						if (!occurMap.containsKey(t[4])) occurMap.put(t[4], new HashSet<>());
						occurMap.get(t[4]).add(ng.nameFromPath(dir.getPath()));
						direcMap.put(t[4], (int) Math.signum(Double.valueOf(t[5])));
						if (!occurMap.containsKey(t[7])) occurMap.put(t[7], new HashSet<>());
						occurMap.get(t[7]).add(ng.nameFromPath(dir.getPath()));
						direcMap.put(t[7], (int) Math.signum(Double.valueOf(t[8])));
					}
				});
			}
		});

		BufferedWriter sifWriter = FileUtil.newBufferedWriter(outSIFNoExt + ".sif");
		rows.values().stream().map(t -> t[0] + "\t" + t[1] + "\t" + t[2] + "\t\t" + t[3]).distinct().forEach(l -> FileUtil.writeln(l, sifWriter));
		sifWriter.close();

		Integer max = occurMap.values().stream().map(Set::size).max(Integer::compareTo).get();

		ValToColor vtc = new ValToColor(new double[]{-max, 0 , max}, new Color[]{new Color(50, 100, 200), Color.WHITE, new Color(200, 80, 50)});
		SiteEffectCollective sec = new SiteEffectCollective();

		BufferedWriter fmtWriter = FileUtil.newBufferedWriter(outSIFNoExt + ".format");
		occurMap.forEach((id, cases) ->
		{
			String mod = getModLetter(id);
			Integer eff = sec.getEffectFromID(id);

			if (mod != null) FileUtil.writeln("node\t" + id.substring(0, id.indexOf("-")) + "\trppasite\t" + id + "|" +
				mod + "|" + vtc.getColorInString(occurMap.get(id).size() * direcMap.get(id)) + "|" +
				(eff == null || eff == 0 ? "0 0 0" : eff == 1 ? "0 180 20" : "180 0 20") + "|" +
				occurMap.get(id).stream().sorted().collect(Collectors.toList()), fmtWriter);
			else
			{
				FileUtil.writeln("node\t" + id + "\tcolor\t" + vtc.getColorInString(occurMap.get(id).size() * direcMap.get(id)), fmtWriter);
				FileUtil.writeln("node\t" + id + "\ttooltip\t" + occurMap.get(id).stream().sorted().collect(Collectors.toList()), fmtWriter);
			}
		});

		fmtWriter.close();
	}

	private static String getModLetter(String id)
	{
		String[] t = id.split("-");
		for (String s : t)
		{
			if (s.equals("P")) return "p";
			if (s.equals("A")) return "a";
			if (s.contains("active")) return null;
//			if (s.equals("active")) return "!";
//			if (s.equals("inactive")) return "i";
		}
		return null;
	}

	private static Map<String, Integer> getRBFocus()
	{
		Map<String, Integer> map = new HashMap<>();
		map.put("RB1-S795-S807-P", 1);
		map.put("RB1-S795-P-2", 1);
		map.put("RB1-S807-S811-P", 1);
		map.put("RB1-S249-P", 1);
		map.put("RB1-S780-P", 1);
		map.put("RB1-S249-T252-P", 1);
		map.put("RB1-S795-P", 1);
		map.put("RB1-S807-P", 1);
//		map.put("RB1-T821-P", -1);
		map.put("RB1-S780-S788-P", 1);
//		map.put("RB1-S608-S612-P", -1);
		map.put("RB1-T373-P", 1);
		map.put("RB1-T826-P", 1);
		map.put("RB1-Y805-S807-P", 1);
//		map.put("RB1-S612-P", -1);
		map.put("RB1-S780-S788-P-2", 1);
		return map;
	}

	private static Map<String, Integer> getGRKFocus()
	{
		Map<String, Integer> map = new HashMap<>();
		map.put("GRK2-S670-P", -1);
		return map;
	}



	interface NameGetter
	{
		String nameFromPath(String path);
	}

}
