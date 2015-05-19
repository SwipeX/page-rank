import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * Created by tim on 5/18/15.
 */
public class Strawman {
    public static void main(String[] args) throws IOException, InterruptedException {
        for (int i = 2; i < 97; i += 2) {
            Runtime rt = Runtime.getRuntime();
            Process proc = rt.exec(new String[]{"srun", "-N2", "-n", "" + i, "./pagerank-serial 1000000"});

            BufferedReader stdInput = new BufferedReader(new
                    InputStreamReader(proc.getInputStream()));

            BufferedReader stdError = new BufferedReader(new
                    InputStreamReader(proc.getErrorStream()));

// read the output from the command
            System.out.println("Here is the standard output of the command:\n");
            String s = null;
            while ((s = stdInput.readLine()) != null) {
                System.out.println(s);
            }

// read any errors from the attempted command
            System.out.println("Here is the standard error of the command (if any):\n");
            while ((s = stdError.readLine()) != null) {
                System.out.println(s);
            }
        }
    }
}
