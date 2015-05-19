import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * Created by tim on 5/18/15.
 */
public class Strawman {
    public static void main(String[] args) throws IOException {
        Runtime rt = Runtime.getRuntime();
        for (int i = 2; i < 97; i++) {
            Process proc = rt.exec(String.format("srun -N2 -n %d ./pagerank-serial 1000000", i));

            BufferedReader stdInput = new BufferedReader(new
                    InputStreamReader(proc.getInputStream()));

            String read;
            while ((read = stdInput.readLine()) != null) {
                System.out.println("#" + read);
            }
        }
    }
}
