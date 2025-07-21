// Navigation functionality for Spacemake reports

document.addEventListener("DOMContentLoaded", function() {
    generateNavigation();
    setupScrollSpy();
});

function generateNavigation() {
    const navItems = document.getElementById("navItems");
    if (!navItems) return;
    
    navItems.innerHTML = "";
    
    // Find all headers in the notebook
    const headers = document.querySelectorAll(".jp-RenderedMarkdown h1, .jp-RenderedMarkdown h2, .jp-RenderedMarkdown h3");
    
    headers.forEach((header, index) => {
        const id = "nav-header-" + index;
        header.id = id;
        
        const navItem = document.createElement("div");
        navItem.className = "nav-item";
        
        // Add level class for indentation
        if (header.tagName === "H2") {
            navItem.classList.add("level-2");
        } else if (header.tagName === "H3") {
            navItem.classList.add("level-3");
        }
        
        // Clean up header text
        let headerText = header.textContent.replace(/^#+\s*/, "").trim();
        if (headerText.length > 50) {
            headerText = headerText.substring(0, 47) + "...";
        }
        navItem.textContent = headerText;
        
        navItem.onclick = function() {
            // Remove active class from all nav items
            document.querySelectorAll(".nav-item").forEach(item => {
                item.classList.remove("active");
            });
            
            // Add active class to clicked item
            this.classList.add("active");
            
            // Scroll to header with offset for fixed elements
            const headerRect = header.getBoundingClientRect();
            const absoluteElementTop = headerRect.top + window.pageYOffset;
            const offset = 20; // Offset from top
            
            window.scrollTo({
                top: absoluteElementTop - offset,
                behavior: "smooth"
            });
            
            // Close mobile sidebar if open
            if (window.innerWidth <= 768) {
                document.getElementById("navSidebar").classList.remove("show");
            }
        };
        
        navItems.appendChild(navItem);
    });
}

function setupScrollSpy() {
    let ticking = false;
    
    function updateActiveNavItem() {
        const headers = document.querySelectorAll(".text_cell_render h1, .text_cell_render h2, .text_cell_render h3");
        const navItems = document.querySelectorAll(".nav-item");
        
        let current = "";
        let currentIndex = -1;
        
        headers.forEach((header, index) => {
            const rect = header.getBoundingClientRect();
            if (rect.top <= 100) {
                current = header.id;
                currentIndex = index;
            }
        });
        
        // Update active nav item
        navItems.forEach(item => item.classList.remove("active"));
        if (currentIndex >= 0 && navItems[currentIndex]) {
            navItems[currentIndex].classList.add("active");
        }
        
        ticking = false;
    }
    
    function requestTick() {
        if (!ticking) {
            requestAnimationFrame(updateActiveNavItem);
            ticking = true;
        }
    }
    
    window.addEventListener("scroll", requestTick);
}

function toggleSidebar() {
    const sidebar = document.getElementById("navSidebar");
    sidebar.classList.toggle("show");
}

// Close sidebar when clicking outside on mobile
document.addEventListener("click", function(event) {
    const sidebar = document.getElementById("navSidebar");
    const toggle = document.querySelector(".mobile-nav-toggle");
    
    if (window.innerWidth <= 768 && 
        !sidebar.contains(event.target) && 
        !toggle.contains(event.target)) {
        sidebar.classList.remove("show");
    }
});