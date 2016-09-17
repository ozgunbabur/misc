package org.panda.misc;

import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Stream;

/**
 * Created by babur on 4/22/16.
 */
public class MutexReader
{
	public static Set<Group> readMutexResults(String dir)
	{
		String file = dir + "/ranked-groups.txt";
		return readResults(dir, file);
	}

	public static Set<Group> readCoocResults(String dir)
	{
		String file = dir + "/cooc-groups.txt";
		return readResults(dir, file);
	}

	private static Set<Group> readResults(String dir, String file)
	{
		boolean hasQval = hasQVal(file);

		Set<Group> groups = new HashSet<>();
		readGroups(file, groups, hasQval, dir);

		return groups;
	}

	public static Map<String, Double> readBestScores(String dir)
	{
		return convertGroupsToGeneBestScores(readMutexResults(dir));
	}

	public static Map<String, Double> readBestScoresRecursive(String dir)
	{
		return convertGroupsToGeneBestScores(readMutexResultsRecursive(dir, new HashSet<>()));
	}

	public static Map<String, Double> convertGroupsToGeneBestScores(Set<Group> groups)
	{
		Map<String, Double> scores = new HashMap<>();
		for (Group group : groups)
		{
			for (String gene : group.genes)
			{
				if (!scores.containsKey(gene) || scores.get(gene) > group.score) scores.put(gene, group.score);
			}
		}
		return scores;
	}

	/**
	 * Pairs are keys of the result map, where key is "g1 g2" and g1.compareTo(g2) < 0.
	 */
	public static Map<String, Double> convertGroupsToPairBestScores(Set<Group> groups)
	{
		Map<String, Double> scores = new HashMap<>();
		for (Group group : groups)
		{
			for (String g1 : group.genes)
			{
				for (String g2 : group.genes)
				{
					if (g1.compareTo(g2) < 0)
					{
						String key = g1 + " " + g2;
						if (!scores.containsKey(key) || scores.get(key) > group.score) scores.put(key, group.score);
					}
				}
			}
		}
		return scores;
	}

	public static boolean hasQVal(String path)
	{
		try(Stream<String> stream = Files.lines(Paths.get(path)))
		{
			return stream.limit(1).anyMatch(line -> line.split("\t").length >= 3);
		}
		catch (IOException e)
		{
			throw new AssertionError(e);
		}
	}

	public static void readGroups(String path, Set<Group> groups, boolean hasQval, String dir)
	{
		try(Stream<String> stream = Files.lines(Paths.get(path)))
		{
			stream.skip(1).filter(l -> !l.isEmpty()).map(line -> line.split("\t")).forEach(token ->
				groups.add(new Group(Arrays.asList(token).subList(hasQval ? 2 : 1, token.length),
					dir, Double.parseDouble(token[0]))));
		}
		catch (IOException e)
		{
			throw new AssertionError(e);
		}
	}

	public static Set<Group> readMutexResultsRecursive(String dir, Set<Group> result)
	{
		if (new File(dir + "/ranked-groups.txt").exists())
			result.addAll(readMutexResults(dir));

		Arrays.stream(new File(dir).listFiles()).filter(File::isDirectory).forEach(d ->
			readMutexResultsRecursive(d.getPath(), result));

		return result;
	}

	public static Set<Group> readMutexResultsRecursive(String dir, Set<Group> result, DirectoryFilter filter)
	{
		if (new File(dir + "/ranked-groups.txt").exists())
		{
			result.addAll(readMutexResults(dir));
		}

		Arrays.stream(new File(dir).listFiles()).filter(File::isDirectory).filter(filter::process).forEach(d ->
			readMutexResultsRecursive(d.getPath(), result, filter));

		return result;
	}

	public static class Group
	{
		public List<String> genes;
		public String fromDir;
		public double score;

		public Group(List<String> genes, String fromDir, double score)
		{
			this.genes = genes;
			this.fromDir = fromDir;
			this.score = score;
		}

		public double getScore()
		{
			return score;
		}

		@Override
		public String toString()
		{
			return score + "\t" + CollectionUtil.merge(genes, " ") + "\t" + fromDir;
		}

		public void shortenLoc(String remove)
		{
			if (fromDir.startsWith(remove)) fromDir = fromDir.substring(remove.length());
		}

		public String geneKey()
		{
			StringBuilder sb = new StringBuilder();
			genes.stream().sorted().forEach(g -> sb.append("\t").append(g));
			return sb.toString().trim();
		}

		@Override
		public int hashCode()
		{
			int h = 0;
			for (String gene : genes)
			{
				h += gene.hashCode();
			}
			h += fromDir.hashCode();
			h += Double.hashCode(score);
			return h;
		}

		@Override
		public boolean equals(Object obj)
		{
			return obj instanceof Group && ((Group) obj).genes.containsAll(genes) &&
				((Group) obj).genes.size() == genes.size() && ((Group) obj).fromDir.equals(fromDir) &&
				((Group) obj).score == score;
		}
	}

	public interface DirectoryFilter
	{
		boolean process(File f);
	}

	public static void main(String[] args)
	{
		String dir = "/home/babur/Downloads/PanCan/";
		readMutexResultsRecursive(dir, new HashSet<>()).stream()
			.sorted(Comparator.comparing(Group::getScore)).filter(g -> g.score < 0.01).peek(g -> g.shortenLoc(dir))
			.map(g -> g.toString().replaceAll("/", "\t")).forEach(System.out::println);
	}
}
