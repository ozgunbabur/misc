package org.panda.misc.teaching;

import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.TermCounter;
import org.panda.utility.statistics.Summary;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * This is for generating a weighted average grade from the individual grades of each zybook assignment.
 */
public class Grader2
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/Users/ozgun/Documents/Teaching/CS220/grades-spring-2021/";
		Map<String, Double> grades = getOverallGrades(dir + "CS_220_Spring_2021_grades.csv", dir + "late-list.txt", dir + "UMBCS220BaburSpring2021_assignment_report_2021-04-22_1620.csv");

		TermCounter tc = new TermCounter();
		grades.keySet().stream().sorted(Comparator.comparing(o -> o.substring(o.lastIndexOf(" ") + 1))).peek(name -> tc.addTerm(getLetterGrade(grades.get(name)))).forEach(name -> System.out.println(getLetterGrade(grades.get(name)) + "\t" + grades.get(name) + "\t" + name));
		System.out.println();
		tc.print();
	}

	private static Map<String, Double> readZyBookGrades(String inFile, String outFile) throws IOException
	{
		String[] header = Files.lines(Paths.get(inFile)).findFirst().get().split(",");
		double[] weights = new double[header.length];

		for (int i = 3; i < header.length; i++)
		{
			weights[i] = Double.parseDouble(header[i].substring(header[i].indexOf("(") + 1, header[i].indexOf(")")));
		}

		double totW = Summary.sum(weights);
		for (int i = 0; i < weights.length; i++)
		{
			weights[i] /= totW;
		}

		Map<String, Double> overall = new HashMap<>();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write(header[0] + "\t" + header[1] + "\tGrade");

		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split(",")).forEach(t ->
		{
			FileUtil.lnwrite(t[0] + "\t" + t[1] + "\t", writer);

			double tot = 0;
			for (int i = 3; i < weights.length; i++)
			{
				tot += Double.parseDouble(t[i]) * weights[i];
			}

			overall.put(t[1] + " " + t[0], tot);
			FileUtil.write("" + (int) Math.round(tot), writer);
		});

		writer.close();
		return overall;
	}

	private static Map<String, Double> getOverallGrades(String gradescopeFile, String latesFile, String zyBookFile) throws IOException
	{
		Map<String, Double> zyBookGrades = readZyBookGrades(zyBookFile, zyBookFile + "-summed.txt");
//		Map<String, Set<String>> latesMap = getLateAssignments(latesFile);

		String[] header = Files.lines(Paths.get(gradescopeFile)).findFirst().get().split(",");

		Set<String> homeworks = Arrays.stream(header).filter(n -> n.startsWith("Homework") && !n.contains(" - ")).collect(Collectors.toSet());
		Set<String> quizzes = Arrays.stream(header).filter(n -> n.startsWith("End") && n.split(" - ").length < 3).collect(Collectors.toSet());

		Map<String, Double> weightMap = new HashMap<>();
		homeworks.forEach(n -> weightMap.put(n, 0.15 / homeworks.size()));
		quizzes.forEach(n -> weightMap.put(n, 1D / quizzes.size())); // Since quizzes are graded over 1

		weightMap.put("Midterm Exam", 0.3);
		weightMap.put("Final Exam", 0.4);

		Map<String, Double> curveMap = new HashMap<>();
//		curveMap.put("Midterm Exam", 9D);
//		curveMap.put("Final Exam", 15D);

//		Map<String, Double> maxMap = new HashMap<>();

		Map<String, Double> gradesMap = new HashMap<>();

		Files.lines(Paths.get(gradescopeFile)).skip(1).map(l -> l.split(",")).forEach(t ->
		{
			String name = t[0] + " " + t[1];
			if (name.equals("Hypothetical Student")) return;

			if (!zyBookGrades.containsKey(name))
			{
				System.out.println("not in zybook = " + name);
			}

			double grade = 0;

			double bookGrade = zyBookGrades.getOrDefault(name, 0D) * 0.14;
			grade += bookGrade;

			double quizGrade = 0;
			for (String quiz : quizzes)
			{
				quizGrade += getAssignmentGrade(quiz, t, header, weightMap, curveMap);
			}

			grade += quizGrade;

			double hwGrade = 0;
			for (String homework : homeworks)
			{
				hwGrade += getAssignmentGrade(homework, t, header, weightMap, curveMap);
			}

			grade += hwGrade;

			grade += getAssignmentGrade("Midterm Exam", t, header, weightMap, curveMap);

			grade += getAssignmentGrade("Final Exam", t, header, weightMap, curveMap);

//			for (String ass : weightMap.keySet())
//			{
//				int index = ArrayUtil.indexOf(header, ass);
//				double x = index >= t.length || t[index].isEmpty() ? 0 : Double.parseDouble(t[index]);
//
////				if (maxMap.getOrDefault(ass, 0D) < x) maxMap.put(ass, x);
//
//				if (curveMap.containsKey(ass)) x += curveMap.get(ass);
//
////				if (ass.equals("Final Exam") && x < 40) System.out.println("failed in final = " + name);
//
////				if (latesMap.containsKey(ass) && latesMap.get(ass).contains(name)) x *= 0.8;
//				x *= weightMap.get(ass);
//				grade += x;
//			}

			gradesMap.put(name, grade);
		});

//		for (String ass : maxMap.keySet())
//		{
//			System.out.println(ass + "\t" + maxMap.get(ass));
//		}

		Set<String> students = gradesMap.keySet();
//		students.forEach(s -> gradesMap.put(s, gradesMap.get(s) * 100D / 60D));
		Double max = students.stream().map(gradesMap::get).filter(g -> g != 100).max(Double::compareTo).get();
		double add = 100 - max;
		System.out.println("add = " + add);
		students.forEach(s -> gradesMap.put(s, gradesMap.get(s) + add));

		return gradesMap;
	}

	private static double getAssignmentGrade(String assName, String[] t, String[] header, Map<String, Double> weightMap, Map<String, Double> curveMap)
	{
		int index = ArrayUtil.indexOf(header, assName);
		double x = index >= t.length || t[index].isEmpty() ? 0 : Double.parseDouble(t[index]);

		if (curveMap.containsKey(assName)) x += curveMap.get(assName);

		x *= weightMap.get(assName);
		return x;
	}

	private static Map<String, Set<String>> getLateAssignments(String file) throws FileNotFoundException
	{
		Map<String, Set<String>> map = new HashMap<>();

		Set<String> names = new HashSet<>();

		Scanner scanner = new Scanner(new File(file));
		while(scanner.hasNextLine())
		{
			String line = scanner.nextLine();
			if (line.isEmpty()) continue;

			if (line.startsWith(">"))
			{
				line = line.substring(1);
				names = new HashSet<>();
				map.put(line, names);
			}
			else names.add(line);
		}

		return map;
	}

//	private static String getLetterGrade(double v)
//	{
//		if (v < 12) return "FAI";
//		else if (v < 21) return "CAU";
//		else return "SAT";
//	}

	private static String getLetterGrade(double v)
	{
		if (v < 40) return "F";
		else if (v < 45) return "D-";
		else if (v < 50) return "D";
		else if (v < 55) return "D+";
		else if (v < 60) return "C-";
		else if (v < 65) return "C";
		else if (v < 70) return "C+";
		else if (v < 75) return "B-";
		else if (v < 80) return "B";
		else if (v < 85) return "B+";
		else if (v < 90) return "A-";
		else return "A";
	}
}
