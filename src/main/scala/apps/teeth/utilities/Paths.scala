package apps.teeth.utilities

import java.io.File

object Paths {
  val localPath = "data/"
  val serverPath = "/export/skulls/projects/teeth/data/reference/LowerMolar/"
  val generalPath = new File(serverPath +"lowermolar1") // Update to SMIR folder

}
